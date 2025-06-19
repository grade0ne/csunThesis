import numpy as np
import matplotlib.pyplot as plt
from skimage.io import imread
from skimage.morphology import skeletonize, binary_closing, remove_small_objects, disk
from skimage.draw import line, polygon2mask
from skimage.transform import resize
from scipy.ndimage import binary_fill_holes, convolve, map_coordinates
from scipy.interpolate import splprep, splev
import networkx as nx
import plotly.graph_objects as go

# ========== Preprocessing Functions ==========

def upscale_binary_mask(mask, scale=2, threshold=0.5):
    upscaled = resize(mask.astype(float),
                      (mask.shape[0] * scale, mask.shape[1] * scale),
                      order=1, preserve_range=True, anti_aliasing=True)
    return upscaled > threshold

def clean_binary_mask(mask, close_footprint_size = 13, rm_objects_size = 100):
    cleaned = binary_closing(mask, footprint = disk(close_footprint_size))
    cleaned = remove_small_objects(cleaned, min_size = rm_objects_size)
    cleaned = binary_fill_holes(cleaned)
    return cleaned

# ========== Skeleton & Midline Extraction ==========

def extract_midline_path(skeleton):
    neighbor_kernel = np.array([[1, 1, 1],
                                [1, 0, 1],
                                [1, 1, 1]])

    neighbor_count = convolve(skeleton.astype(np.uint8),
                              neighbor_kernel,
                              mode='constant',
                              cval=0)
    endpoints = np.logical_and(skeleton, neighbor_count == 1)
    endpoint_coords = np.column_stack(np.where(endpoints))
    if len(endpoint_coords) < 2:
        raise ValueError("Midline has fewer than two endpoints.")

    G = nx.Graph()
    ys, xs = np.where(skeleton)
    for y, x in zip(ys, xs):
        for dy in [-1, 0, 1]:
            for dx in [-1, 0, 1]:
                if dy == 0 and dx == 0:
                    continue
                yn, xn = y + dy, x + dx
                if 0 <= yn < skeleton.shape[0] and 0 <= xn < skeleton.shape[1]:
                    if skeleton[yn, xn]:
                        G.add_edge((y, x), (yn, xn))

    start = tuple(endpoint_coords[0])
    end = tuple(endpoint_coords[1])
    return nx.shortest_path(G, source=start, target=end)

# ========== Spline and Midline Utilities ==========

def fit_midline_spline(path, smoothing=0, spacing=3.0):
    path = np.array(path)
    x, y = path[:, 1], path[:, 0]
    dx, dy = np.diff(x), np.diff(y)
    distances = np.sqrt(dx**2 + dy**2)
    arc_length = np.concatenate([[0], np.cumsum(distances)])
    arc_length_norm = arc_length / arc_length[-1]
    tck, _ = splprep([x, y], u=arc_length_norm, s=smoothing)

    u_dense = np.linspace(0, 1, 1000)
    x_dense, y_dense = splev(u_dense, tck)
    points = np.column_stack((x_dense, y_dense))
    step_lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    arc_lengths = np.concatenate([[0], np.cumsum(step_lengths)])
    total_length = arc_lengths[-1]

    num_samples = max(int(total_length / spacing), 2)
    target_lengths = np.linspace(0, total_length, num_samples)
    x_resampled = np.interp(target_lengths, arc_lengths, x_dense)
    y_resampled = np.interp(target_lengths, arc_lengths, y_dense)

    return list(zip(np.round(y_resampled).astype(int), np.round(x_resampled).astype(int)))

def deduplicate_rounded_path(path):
    seen = set()
    cleaned = []
    for pt in path:
        rounded = tuple(np.round(pt).astype(int))
        if rounded not in seen:
            seen.add(rounded)
            cleaned.append(pt)
    return cleaned


def extend_spline_to_mask(mask, spline_path, spacing=3.0, smoothing_window=5):
    
    def get_smoothed_direction(points, from_start=True):
        points = np.array(points)
        if from_start:
            segment = points[:smoothing_window]
            deltas = segment[1:] - segment[:-1]
        else:
            segment = points[-smoothing_window:]
            deltas = segment[1:] - segment[:-1]
        direction = np.mean(deltas, axis=0)
        norm = np.linalg.norm(direction)
        return direction / norm if norm > 0 else np.array([0.0, 0.0])

    def extend_to_mask_edge(mask, point, direction, max_steps=20):
        y0, x0 = point
        dy, dx = direction
        h, w = mask.shape
        for step in range(1, max_steps + 1):
            y = int(round(y0 + dy * step))
            x = int(round(x0 + dx * step))
            if not (0 <= y < h and 0 <= x < w):
                break
            if mask[y, x] == 0:
                y = int(round(y0 + dy * (step - 1)))
                x = int(round(x0 + dx * (step - 1)))
                return (y, x)
        return None

    def interpolate_extension(p1, p2, spacing=3.0):
        p1 = np.array(p1, dtype=float)
        p2 = np.array(p2, dtype=float)
        extension_length = np.linalg.norm(p2 - p1)
        num_points = max(int(np.round(extension_length / spacing)), 1)
        return [tuple(p1 + (p2 - p1) * t)
                for t in np.linspace(0, 1, num_points + 1)[1:]]

    padded_spline = []

    start_dir = get_smoothed_direction(spline_path, from_start=True)
    end_dir = get_smoothed_direction(spline_path, from_start=False)

    start_pad = extend_to_mask_edge(mask, spline_path[0], start_dir, max_steps=20)
    end_pad = extend_to_mask_edge(mask, spline_path[-1], end_dir, max_steps=20)

    if start_pad:
        start_interp = interpolate_extension(spline_path[0], start_pad, spacing=spacing)
        padded_spline = start_interp[::-1] + padded_spline  # prepend reversed

    padded_spline += spline_path

    if end_pad:
        end_interp = interpolate_extension(spline_path[-1], end_pad, spacing=spacing)
        padded_spline.extend(end_interp)

    return deduplicate_rounded_path(padded_spline)


# ========== Geometry and Sampling ==========

def get_unit_direction(p1, p2):
    dy, dx = p2[0] - p1[0], p2[1] - p1[1]
    norm = np.hypot(dy, dx)
    return (dy / norm, dx / norm) if norm != 0 else (0, 0)

def compute_normals(path, offset=2):
    path = np.array(path).astype(float)
    tangents = np.zeros_like(path)
    for i in range(len(path)):
        i1 = max(0, i - offset)
        i2 = min(len(path) - 1, i + offset)
        delta = path[i2] - path[i1]
        norm = np.linalg.norm(delta)
        tangents[i] = delta / norm if norm > 1e-6 else [0, 0]
    normals = np.zeros_like(tangents)
    normals[:, 0] = -tangents[:, 1]
    normals[:, 1] = tangents[:, 0]
    return normals

def compute_arc_lengths(path):
    path = np.array(path)
    diffs = np.diff(path, axis=0)
    return np.concatenate([[0], np.cumsum(np.linalg.norm(diffs, axis=1))])

# ========== Width Sampling ==========

def sample_widths(mask, midline, normals, max_distance=40, step_size=0.5):
    mask = mask.astype(np.uint8)
    widths, edge_points = [], []
    for (y0, x0), (dy, dx) in zip(midline, normals):
        norm = np.hypot(dy, dx)
        if norm == 0:
            widths.append(np.nan)
            edge_points.append((None, None))
            continue
        dy, dx = dy / norm, dx / norm
        found1 = found2 = False
        for d in np.arange(0, max_distance, step_size):
            y_pos, x_pos = y0 + dy * d, x0 + dx * d
            y_neg, x_neg = y0 - dy * d, x0 - dx * d
            val_pos = map_coordinates(mask, [[y_pos], [x_pos]], order=1, mode='constant')[0]
            val_neg = map_coordinates(mask, [[y_neg], [x_neg]], order=1, mode='constant')[0]
            if not found1 and val_pos < 0.5:
                p1 = (y0 + dy * (d - step_size), x0 + dx * (d - step_size))
                found1 = True
            if not found2 and val_neg < 0.5:
                p2 = (y0 - dy * (d - step_size), x0 - dx * (d - step_size))
                found2 = True
            if found1 and found2:
                break
        if found1 and found2:
            dist = np.hypot(p1[0] - p2[0], p1[1] - p2[1])
            widths.append(dist)
            edge_points.append((p1, p2))
        else:
            widths.append(np.nan)
            edge_points.append((None, None))
    return widths, edge_points

# ========== 3D Ring Visualization ==========

def plot_oriented_3d_wireframe(spline_coords, arc_lengths, widths, step=2, n_points=40):
    fig = go.Figure()
    spline_coords = np.array(spline_coords).astype(float)
    arc_lengths = np.array(arc_lengths)
    centers = np.column_stack((spline_coords[:, 1], spline_coords[:, 0], arc_lengths))

    tangents = np.gradient(centers, axis=0)
    tangents /= np.linalg.norm(tangents, axis=1, keepdims=True)

    for i in range(0, len(widths), step):
        if np.isnan(widths[i]):
            continue
        T = tangents[i]
        center = centers[i]
        ref = np.array([0, 0, 1]) if abs(T[2]) < 0.99 else np.array([1, 0, 0])
        N = np.cross(T, ref)
        N /= np.linalg.norm(N)
        B = np.cross(T, N)
        B /= np.linalg.norm(B)
        r = widths[i] / 2
        theta = np.linspace(0, 2 * np.pi, n_points)
        ring = np.outer(np.cos(theta), N) + np.outer(np.sin(theta), B)
        ring = center + r * ring
        fig.add_trace(go.Scatter3d(
            x=ring[:, 0], y=ring[:, 1], z=ring[:, 2],
            mode='lines',
            line=dict(color='royalblue', width=2),
            showlegend=False
        ))

    fig.update_layout(
        scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Arc Length (Z)', aspectmode='data'),
        title='3D Rotifer Reconstruction (Stable Frame)',
        margin=dict(l=0, r=0, t=40, b=0)
    )
    fig.show()

# ========== Visualization ========== 

def plot_mask_comparison(original_mask, cleaned_mask, skeleton=None):
    
    if original_mask.shape != cleaned_mask.shape:
        raise ValueError("Original and cleaned masks must have the same shape for comparison visualization.")

    overlay = np.zeros((*original_mask.shape, 3), dtype=np.uint8)
    overlay[..., 0] = ((cleaned_mask == 1) & (original_mask == 0)) * 255  # red: added fill
    overlay[..., 1] = ((original_mask == 1) & (cleaned_mask == 1)) * 255  # green: preserved
    overlay[..., 2] = ((original_mask == 1) & (cleaned_mask == 0)) * 255  # blue: removed

    if skeleton is not None:
        yx = np.argwhere(skeleton)
        for y, x in yx:
            overlay[y, x] = [255, 255, 0]  # yellow overlay for skeleton

    plt.figure(figsize=(6, 6))
    plt.imshow(overlay)
    plt.title("Mask Comparison (Red: Fill-ins, Blue: Spills, Yellow: Skeleton)")
    plt.axis('off')
    plt.tight_layout()
    plt.show()

def plot_mask_with_normals(mask, spline, normals, step=1):
    
    plt.figure(figsize=(6, 6))
    plt.imshow(mask, cmap='gray')
    yy, xx = zip(*spline)
    plt.plot(xx, yy, 'r-', label='Spline')

    for (y, x), (dy, dx) in zip(spline[::step], normals[::step]):
        plt.arrow(x, y, dx*5, dy*5, color='cyan', head_width=1)

    plt.title("Spline and Normals")
    plt.axis('off')
    plt.tight_layout()
    plt.show()

def plot_3d_rings_with_estimates(midline, widths, arcs, step=2, n_points=40):
    """
    Plot 3D rings for measured widths and overlay frustum and r-avg estimate rings,
    correctly oriented orthogonal to the spline using robust local frames.
    Also includes a 3D spline vector trace.
    """
    fig = go.Figure()
    midline = np.array(midline).astype(float)
    centers = np.column_stack((midline[:, 1], midline[:, 0], arcs))  # (x, y, z)

    # Plot the midline spline as a 3D vector trace
    fig.add_trace(go.Scatter3d(
        x=centers[:, 0], y=centers[:, 1], z=centers[:, 2],
        mode='lines',
        line=dict(color='black', width=4),
        name='Midline Spline'
    ))

    # Compute local frames using orthogonal basis from cross products
    tangents = np.gradient(centers, axis=0)
    tangents /= np.linalg.norm(tangents, axis=1, keepdims=True)
    frames = []

    for i, T in enumerate(tangents):
        T = T / np.linalg.norm(T)
        ref = np.array([0, 0, 1]) if abs(T[2]) < 0.99 else np.array([1, 0, 0])
        N = np.cross(T, ref)
        N /= np.linalg.norm(N)
        B = np.cross(T, N)
        B /= np.linalg.norm(B)
        frames.append((T, N, B))

    # Draw rings
    for i in range(0, len(widths) - 1, step):
        if np.isnan(widths[i]) or np.isnan(widths[i + 1]):
            continue
        center = centers[i]
        T, N, B = frames[i]

        # Get true radii
        r_measured = widths[i] / 2
        r_frustum = (widths[i] + widths[i + 1]) / 4
        r_ravg = np.sqrt((widths[i] / 2)**2 + (widths[i + 1] / 2)**2) / 2

        for radius, color, label in zip([r_measured, r_frustum, r_ravg],
                                        ['royalblue', 'red', 'green'],
                                        ['Measured', 'Frustum', 'r_avg']):
            theta = np.linspace(0, 2 * np.pi, n_points)
            x_ring = center[0] + radius * (np.cos(theta) * N[0] + np.sin(theta) * B[0])
            y_ring = center[1] + radius * (np.cos(theta) * N[1] + np.sin(theta) * B[1])
            z_ring = center[2] + radius * (np.cos(theta) * N[2] + np.sin(theta) * B[2])

            fig.add_trace(go.Scatter3d(
                x=x_ring, y=y_ring, z=z_ring,
                mode='lines',
                line=dict(color=color, width=2),
                name=label if i == 0 else None,
                showlegend=(i == 0)
            ))

    fig.update_layout(
        scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Arc Length', aspectmode='data'),
        title='3D Width Rings with Volume Estimations (Properly Oriented)',
        margin=dict(l=0, r=0, t=40, b=0)
    )
    fig.show()

def plot_3d_volume_mesh(midline, widths, arcs, step=2, n_points=40):
    """
    Generate a 3D volume mesh surface by connecting cross-sectional rings along the midline.
    """
    import plotly.graph_objects as go
    midline = np.array(midline).astype(float)
    centers = np.column_stack((midline[:, 1], midline[:, 0], arcs))

    # Compute frames
    tangents = np.gradient(centers, axis=0)
    tangents /= np.linalg.norm(tangents, axis=1, keepdims=True)
    frames = []
    for i, T in enumerate(tangents):
        T = T / np.linalg.norm(T)
        ref = np.array([0, 0, 1]) if abs(T[2]) < 0.99 else np.array([1, 0, 0])
        N = np.cross(T, ref)
        N /= np.linalg.norm(N)
        B = np.cross(T, N)
        B /= np.linalg.norm(B)
        frames.append((T, N, B))

    # Build ring points
    ring_list = []
    theta = np.linspace(0, 2 * np.pi, n_points, endpoint=False)
    for i in range(0, len(widths), step):
        if np.isnan(widths[i]):
            continue
        center = centers[i]
        _, N, B = frames[i]
        radius = widths[i] / 2
        ring = [center + radius * (np.cos(t) * N + np.sin(t) * B) for t in theta]
        ring_list.append(ring)

    # Flatten vertices
    vertices = np.vstack(ring_list)
    x, y, z = vertices[:, 0], vertices[:, 1], vertices[:, 2]

    # Create faces between adjacent rings
    faces = []
    for i in range(len(ring_list) - 1):
        for j in range(n_points):
            next_j = (j + 1) % n_points
            a = i * n_points + j
            b = i * n_points + next_j
            c = (i + 1) * n_points + next_j
            d = (i + 1) * n_points + j
            faces.append([a, b, c])
            faces.append([a, c, d])

    i, j, k = zip(*faces)

    fig = go.Figure(data=[go.Mesh3d(
        x=x, y=y, z=z,
        i=i, j=j, k=k,
        color='lightblue', opacity=0.5,
        name='Volume Mesh'
    )])

    fig.update_layout(
        scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Arc Length', aspectmode='data'),
        title='3D Volume Mesh from Width Samples',
        margin=dict(l=0, r=0, t=40, b=0)
    )
    fig.show()

def plot_3d_frustum_mesh(midline, widths, arcs, n_points=40, max_rings=200):
    """
    Generate a 3D mesh with true frustum (truncated cone) segments between each pair of slices.
    Downsamples if necessary to avoid overly large meshes.
    """
    import plotly.graph_objects as go
    midline = np.array(midline).astype(float)
    widths = np.array(widths)
    arcs = np.array(arcs)

    # Filter valid indices
    valid = [i for i in range(len(widths)) if not np.isnan(widths[i])]
    if len(valid) < 2:
        print("Not enough valid width samples to build frustum mesh.")
        return

    # Downsample if necessary
    step = max(1, len(valid) // max_rings)
    indices = valid[::step]
    if len(indices) < 2:
        print("Too few points after downsampling.")
        return

    centers = np.column_stack((midline[indices, 1], midline[indices, 0], arcs[indices]))
    radii = widths[indices] / 2

    # Build frames
    tangents = np.gradient(centers, axis=0)
    tangents /= np.linalg.norm(tangents, axis=1, keepdims=True)
    frames = []
    for i, T in enumerate(tangents):
        ref = np.array([0, 0, 1]) if abs(T[2]) < 0.99 else np.array([1, 0, 0])
        N = np.cross(T, ref)
        N /= np.linalg.norm(N)
        B = np.cross(T, N)
        B /= np.linalg.norm(B)
        frames.append((T, N, B))

    # Build rings
    ring_list = []
    theta = np.linspace(0, 2 * np.pi, n_points, endpoint=False)
    for center, r, (_, N, B) in zip(centers, radii, frames):
        ring = [center + r * (np.cos(t) * N + np.sin(t) * B) for t in theta]
        ring_list.append(ring)

    vertices = np.vstack(ring_list)
    x, y, z = vertices[:, 0], vertices[:, 1], vertices[:, 2]

    faces = []
    for i in range(len(ring_list) - 1):
        for j in range(n_points):
            nj = (j + 1) % n_points
            a = i * n_points + j
            b = i * n_points + nj
            c = (i + 1) * n_points + nj
            d = (i + 1) * n_points + j
            faces.append([a, b, c])
            faces.append([a, c, d])

    if not faces:
        print("No valid faces to render.")
        return

    i, j, k = zip(*faces)

    fig = go.Figure(data=[go.Mesh3d(
        x=x, y=y, z=z,
        i=i, j=j, k=k,
        color='salmon', opacity=0.6,
        name='Frustum Mesh'
    )])

    fig.update_layout(
        scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Arc Length', aspectmode='data'),
        title='3D Volume Mesh Using Frustum Segments (Capped)',
        margin=dict(l=0, r=0, t=40, b=0)
    )
    fig.show()

def compare_frustum_to_mesh_volume(midline, widths, arcs, n_points=40, max_rings=200):
    """
    Compare the analytical frustum volume to a volumetric estimate from a watertight triangle mesh using trimesh.
    """
    import numpy as np
    import trimesh

    midline = np.array(midline).astype(float)
    widths = np.array(widths)
    arcs = np.array(arcs)

    # Filter valid indices
    valid = [i for i in range(len(widths)) if not np.isnan(widths[i])]
    if len(valid) < 2:
        print("Not enough valid width samples to compare volumes.")
        return

    # Downsample if necessary
    step = max(1, len(valid) // max_rings)
    indices = valid[::step]
    if len(indices) < 2:
        print("Too few points after downsampling.")
        return

    centers = np.column_stack((midline[indices, 1], midline[indices, 0], arcs[indices]))
    radii = widths[indices] / 2

    # Analytical frustum volume
    arc_pos = arcs[indices]
    frustum_volume = 0
    for i in range(len(radii) - 1):
        r1, r2 = radii[i], radii[i + 1]
        h = arc_pos[i + 1] - arc_pos[i]
        frustum_volume += (np.pi * h / 3) * (r1**2 + r1 * r2 + r2**2)

    # Generate mesh vertices and faces
    tangents = np.gradient(centers, axis=0)
    tangents /= np.linalg.norm(tangents, axis=1, keepdims=True)
    frames = []
    for i, T in enumerate(tangents):
        ref = np.array([0, 0, 1]) if abs(T[2]) < 0.99 else np.array([1, 0, 0])
        N = np.cross(T, ref)
        N /= np.linalg.norm(N)
        B = np.cross(T, N)
        B /= np.linalg.norm(B)
        frames.append((T, N, B))

    theta = np.linspace(0, 2 * np.pi, n_points, endpoint=False)
    ring_list = []
    for center, r, (_, N, B) in zip(centers, radii, frames):
        ring = [center + r * (np.cos(t) * N + np.sin(t) * B) for t in theta]
        ring_list.append(ring)

    vertices = np.vstack(ring_list)

    # Build triangle faces
    faces = []
    for i in range(len(ring_list) - 1):
        for j in range(n_points):
            nj = (j + 1) % n_points
            a = i * n_points + j
            b = i * n_points + nj
            c = (i + 1) * n_points + nj
            d = (i + 1) * n_points + j
            faces.append([a, b, c])
            faces.append([a, c, d])

    faces = np.array(faces)
    try:
        mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=True)
        mesh_volume = mesh.volume
    except Exception as e:
        print(f"Trimesh volume calculation failed: {e}")
        mesh_volume = np.nan

    print(f"Frustum formula volume: {frustum_volume:.2f} pixels^3")
    print(f"Mesh-integrated volume (trimesh): {mesh_volume:.2f} pixels^3")
    print(f"Difference: {abs(mesh_volume - frustum_volume):.2f} pixels^3")

def plot_3d_projected_2d_normals(midline, normals, arcs, scale=0.3):
    """
    Visualize 2D image-space normals projected into 3D with Z = arc length.
    Each normal is drawn as a cyan arrow from the midline point.
    """
    import numpy as np
    import plotly.graph_objects as go

    midline = np.array(midline).astype(float)
    normals = np.array(normals).astype(float)
    arcs = np.array(arcs).astype(float)

    # Ensure shapes match
    min_len = min(len(midline), len(normals), len(arcs))
    midline = midline[:min_len]
    normals = normals[:min_len]
    arcs = arcs[:min_len]

    # Build 3D points and projected vectors
    arrow_x, arrow_y, arrow_z = [], [], []
    arrow_u, arrow_v, arrow_w = [], [], []

    for (y, x), (dy, dx), z in zip(midline, normals, arcs):
        arrow_x.append(x)
        arrow_y.append(y)
        arrow_z.append(z)
        arrow_u.append(dx * scale)
        arrow_v.append(dy * scale)
        arrow_w.append(0.0)  # flat in XY plane

    fig = go.Figure()

    fig.add_trace(go.Cone(
        x=arrow_x, y=arrow_y, z=arrow_z,
        u=arrow_u, v=arrow_v, w=arrow_w,
        sizemode="absolute", sizeref=1,
        anchor="tail",
        colorscale=[[0, 'cyan'], [1, 'cyan']],
        showscale=False,
        name="2D Projected Normals"
    ))

    fig.add_trace(go.Scatter3d(
        x=[x for (_, x) in midline],
        y=[y for (y, _) in midline],
        z=arcs,
        mode='lines',
        line=dict(color='black', width=3),
        name='Midline'
    ))

    fig.update_layout(
        scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Arc Length (Z)', aspectmode='data'),
        title='Projected 2D Normals in 3D Space',
        margin=dict(l=0, r=0, t=40, b=0)
    )

    fig.show()

def plot_3d_rings_with_mesh_overlay(midline, widths, arcs, step=2, n_points=40, max_rings=200):
    
    import numpy as np
    import plotly.graph_objects as go

    midline = np.array(midline).astype(float)
    widths = np.array(widths)
    arcs = np.array(arcs)
    centers = np.column_stack((midline[:, 1], midline[:, 0], arcs))

    valid = [i for i in range(len(widths) - 1) if not np.isnan(widths[i]) and not np.isnan(widths[i + 1])]
    if len(valid) < 2:
        print("Not enough valid rings to plot.")
        return

    step_down = max(1, len(valid) // max_rings)
    indices = valid[::step_down]

    fig = go.Figure()

    # Midline trace
    fig.add_trace(go.Scatter3d(
        x=centers[:, 0], y=centers[:, 1], z=centers[:, 2],
        mode='lines',
        line=dict(color='black', width=4),
        name='Midline Spline'
    ))

    ring_list = []
    theta = np.linspace(0, 2 * np.pi, n_points, endpoint=False)

    for i in indices:
        center = centers[i]
        if i == 0:
            tangent = centers[i + 1] - centers[i]
        elif i == len(centers) - 1:
            tangent = centers[i] - centers[i - 1]
        else:
            tangent = centers[i + 1] - centers[i - 1]
        tangent /= np.linalg.norm(tangent)

        # Use global Z as default reference
        ref = np.array([0, 0, 1]) if abs(tangent[2]) < 0.99 else np.array([1, 0, 0])
        N = np.cross(tangent, ref)
        N /= np.linalg.norm(N)
        B = np.cross(tangent, N)
        B /= np.linalg.norm(B)

        r_measured = widths[i] / 2
        r_frustum = (widths[i] + widths[i + 1]) / 4
        r_ravg = np.sqrt((widths[i] / 2)**2 + (widths[i + 1] / 2)**2) / 2

        for radius, color, label in zip([r_measured, r_frustum, r_ravg],
                                        ['royalblue', 'red', 'green'],
                                        ['Measured', 'Frustum', 'r_avg']):
            if np.isnan(radius) or radius <= 0:
                continue
            ring = center + radius * (np.outer(np.cos(theta), N) + np.outer(np.sin(theta), B))
            fig.add_trace(go.Scatter3d(
                x=ring[:, 0], y=ring[:, 1], z=ring[:, 2],
                mode='lines',
                line=dict(color=color, width=2),
                name=label if i == indices[0] else None,
                showlegend=(i == indices[0])
            ))

        # Add for mesh
        ring_list.append(center + r_measured * (np.outer(np.cos(theta), N) + np.outer(np.sin(theta), B)))

    # Mesh
    vertices = np.vstack(ring_list)
    x, y, z = vertices[:, 0], vertices[:, 1], vertices[:, 2]
    faces = []
    for i in range(len(ring_list) - 1):
        for j in range(n_points):
            nj = (j + 1) % n_points
            a = i * n_points + j
            b = i * n_points + nj
            c = (i + 1) * n_points + nj
            d = (i + 1) * n_points + j
            faces.append([a, b, c])
            faces.append([a, c, d])

    if faces:
        i_f, j_f, k_f = zip(*faces)
        fig.add_trace(go.Mesh3d(
            x=x, y=y, z=z,
            i=i_f, j=j_f, k=k_f,
            color='lightblue', opacity=0.5,
            name='Volume Mesh'
        ))

    fig.update_layout(
        scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Arc Length (Z)', aspectmode='data'),
        title='3D Rings + Volume Mesh (Simplified)',
        margin=dict(l=0, r=0, t=40, b=0)
    )

    fig.show()

def plot_voxel_volume_as_cubes_filled(midline, normals, widths, arcs, image_shape, voxel_size=1.0, interp_steps=3):
    """
    Render voxel volume using filled cubic voxels with linear interpolation between slices.
    """
    from skimage.draw import polygon2mask

    midline = np.array(midline)
    normals = np.array(normals)
    widths = np.array(widths)
    arcs = np.array(arcs)

    margin = 10
    ys, xs = midline[:, 0], midline[:, 1]
    min_y, max_y = int(np.floor(ys.min())) - margin, int(np.ceil(ys.max())) + margin
    min_x, max_x = int(np.floor(xs.min())) - margin, int(np.ceil(xs.max())) + margin
    min_z, max_z = 0, int(np.ceil(arcs.max())) + 1

    shape_y = max_y - min_y
    shape_x = max_x - min_x
    shape_z = max_z

    volume = np.zeros((shape_z, shape_y, shape_x), dtype=bool)

    n_points = 40
    theta = np.linspace(0, 2 * np.pi, n_points, endpoint=False)

    for i in range(len(midline) - 1):
        if np.isnan(widths[i]) or np.isnan(widths[i + 1]):
            continue

        c1, c2 = midline[i], midline[i + 1]
        n1, n2 = normals[i], normals[i + 1]
        w1, w2 = widths[i], widths[i + 1]
        z1, z2 = arcs[i], arcs[i + 1]

        for step in range(interp_steps + 1):
            t = step / (interp_steps + 1)
            center = (1 - t) * c1 + t * c2
            normal = (1 - t) * n1 + t * n2
            normal /= np.linalg.norm(normal)
            radius = (1 - t) * w1 + t * w2
            arc = (1 - t) * z1 + t * z2

            B = np.array([-normal[1], normal[0]])
            ring = np.array([
                center + (radius / 2) * (np.cos(a) * normal + np.sin(a) * B)
                for a in theta
            ])

            yy, xx = ring[:, 0] - min_y, ring[:, 1] - min_x
            mask = polygon2mask((shape_y, shape_x), np.column_stack((yy, xx)))

            slice_z = int(np.round(arc))
            if 0 <= slice_z < shape_z:
                volume[slice_z] |= mask

    # Extract voxel positions
    z_idx, y_idx, x_idx = np.nonzero(volume)
    x_coords = x_idx * voxel_size + min_x
    y_coords = y_idx * voxel_size + min_y
    z_coords = z_idx * voxel_size

    # Build mesh
    cube_size = voxel_size
    vertices = []
    faces = []

    for x, y, z in zip(x_coords, y_coords, z_coords):
        base_idx = len(vertices)
        verts = np.array([
            [x, y, z],
            [x + cube_size, y, z],
            [x + cube_size, y + cube_size, z],
            [x, y + cube_size, z],
            [x, y, z + cube_size],
            [x + cube_size, y, z + cube_size],
            [x + cube_size, y + cube_size, z + cube_size],
            [x, y + cube_size, z + cube_size],
        ])
        f = np.array([
            [0, 1, 2], [0, 2, 3],
            [4, 5, 6], [4, 6, 7],
            [0, 1, 5], [0, 5, 4],
            [2, 3, 7], [2, 7, 6],
            [1, 2, 6], [1, 6, 5],
            [0, 3, 7], [0, 7, 4]
        ]) + base_idx
        vertices.extend(verts)
        faces.extend(f)

    vertices = np.array(vertices)
    faces = np.array(faces)
    i, j, k = faces.T

    fig = go.Figure(data=go.Mesh3d(
        x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2],
        i=i, j=j, k=k,
        color='mediumpurple', opacity=0.4,
        showscale=False,
        name='Voxel Mesh (Filled)'
    ))

    fig.update_layout(
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Arc Length (Z)',
            aspectmode='data'
        ),
        title='Filled Voxelized Volume with Interpolated Slices',
        margin=dict(l=0, r=0, t=40, b=0)
    )

    fig.show()


def run_full_visualization(original_mask, cleaned_mask, skeleton, spline, normals, widths, arcs):
    
    # plot_mask_comparison(original_mask, cleaned_mask, skeleton)
    # plot_mask_with_normals(cleaned_mask, spline, normals)
    # plot_3d_rings_with_estimates(spline, widths, arcs)
    # plot_3d_volume_mesh(spline, widths, arcs)
    # plot_3d_frustum_mesh(spline, widths, arcs)
    compare_frustum_to_mesh_volume(spline, widths, arcs)
    plot_3d_rings_with_mesh_overlay(spline, widths, arcs, normals)
    plot_voxel_volume_as_cubes_filled(padded_spline, normals, widths, arcs,
                                        image_shape=cleaned.shape,
                                        voxel_size=1.0,
                                        interp_steps=2)


# ========== Debug Functions ==========

def compare_trimmed_volume_effect(midline, widths, arcs, trim_n=2):
    """
    Compare volume estimates with and without the first `trim_n` samples,
    to assess how much early spline artifacts affect total volume.
    """
    import numpy as np

    def estimate_frustum_volume(radii, arc_positions):
        volume = 0
        for i in range(len(radii) - 1):
            r1, r2 = radii[i], radii[i + 1]
            h = arc_positions[i + 1] - arc_positions[i]
            volume += (np.pi * h / 3) * (r1**2 + r1 * r2 + r2**2)
        return volume

    # Full volume
    valid_full = [i for i in range(len(widths)) if not np.isnan(widths[i])]
    if len(valid_full) < 2:
        print("Not enough valid width samples to compute volume.")
        return

    radii_full = np.array(widths)[valid_full] / 2
    arcs_full = np.array(arcs)[valid_full]
    vol_full = estimate_frustum_volume(radii_full, arcs_full)

    # Trimmed volume
    if len(valid_full) <= trim_n:
        print("Not enough samples after trimming to estimate volume.")
        return

    valid_trimmed = valid_full[trim_n:]
    radii_trimmed = np.array(widths)[valid_trimmed] / 2
    arcs_trimmed = np.array(arcs)[valid_trimmed]
    arcs_trimmed -= arcs_trimmed[0]  # Re-zero arc positions

    vol_trimmed = estimate_frustum_volume(radii_trimmed, arcs_trimmed)

    diff = vol_full - vol_trimmed
    percent = (diff / vol_full) * 100

    print(f"Frustum volume (full):   {vol_full:.2f} pixels³")
    print(f"Frustum volume (trimmed): {vol_trimmed:.2f} pixels³")
    print(f"Difference: {diff:.2f} pixels³ ({percent:.2f}%)")

# ========== Example Usage ==========

if __name__ == "__main__":
    
    # Load and preprocess
    mask_path = "G:/Research/imageAnalysis/partical analysis results/ROImasks/c7-17-1_0031_s_SEG_BINf_roi2_mask.tif"
    original_mask = imread(mask_path) == 0
    upscaled = upscale_binary_mask(mask=original_mask, scale=2)
    cleaned = clean_binary_mask(mask=upscaled, close_footprint_size=13, rm_objects_size=100)

    # Skeleton and midline
    skeleton = skeletonize(cleaned)
    midline_path = extract_midline_path(skeleton)
    spline = fit_midline_spline(midline_path, spacing=3.0)
    padded_spline = extend_spline_to_mask(cleaned, spline, spacing=3.0)

    # Normals and arc lengths
    normals = compute_normals(padded_spline)

    arcs = compute_arc_lengths(padded_spline)

    # Widths
    widths, _ = sample_widths(cleaned, padded_spline, normals)


    compare_trimmed_volume_effect(padded_spline, widths, arcs, trim_n=2)


    # Downscale original mask for comparison if necessary
    if original_mask.shape != cleaned.shape:
        from skimage.transform import resize
        original_mask = resize(original_mask.astype(float), cleaned.shape, order=0, preserve_range=True) > 0.5

    # Run visualizations
    run_full_visualization(original_mask, cleaned, skeleton, padded_spline, normals, widths, arcs)
