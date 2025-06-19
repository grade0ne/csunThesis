import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from skimage.io import imread

from skimage.morphology import skeletonize
from skimage.measure import label, regionprops
from skimage.morphology import binary_closing, binary_opening, remove_small_objects, disk
from skimage.draw import line
from skimage.transform import resize

from scipy.ndimage import binary_fill_holes, convolve, map_coordinates
from scipy.spatial.distance import euclidean
from scipy.interpolate import splprep, splev

import networkx as nx
import plotly.graph_objects as go

import os
from glob import glob

def upscale_binary_mask(mask, scale = 2, threshold = 0.5):
    upscaled = resize(mask.astype(float),
                      (mask.shape[0] * scale, mask.shape[1] * scale),
                      order = 1,
                      preserve_range= True,
                      anti_aliasing= True)
    return upscaled > threshold

example_mask_path = "G:/Research/imageAnalysis/partical analysis results/ROImasks/c7-17-1_0025_s_SEG_BINf_roi0_mask.tif"

binary_img = imread(example_mask_path) == 0      # converts to binary mask where TRUE is rotifer, FALSE is background

binary_up = upscale_binary_mask(binary_img, scale = 2, threshold = 0.5)

cleaned = binary_closing(binary_up, footprint = disk(13))
cleaned = remove_small_objects(cleaned, min_size = 100)
cleaned = binary_fill_holes(cleaned)

skeleton = skeletonize(cleaned)

# Finding endpoints

neighbor_kernel = np.array([[1, 1, 1],
                            [1, 0, 1],
                            [1, 1, 1]])

neighbor_count = convolve(skeleton.astype(np.uint8), 
                          neighbor_kernel, 
                          mode = 'constant', 
                          cval = 0)

endpoints = np.logical_and(skeleton, neighbor_count == 1)
endpoint_coords = np.column_stack(np.where(endpoints))
print(f"Endpoint coords: {endpoint_coords}")

# Midline path

G = nx.Graph()

ys, xs = np.where(skeleton)

for y , x in zip(ys, xs):               # for every ordered pair in the skeleton...
    for dy in [-1, 0, 1]:               # 
        for dx in [-1, 0, 1]:
            if dy == 0 and dx == 0:
                continue
            yn = y + dy
            xn = x + dx

            if 0 <= yn < skeleton.shape[0] and 0 <= xn < skeleton.shape[1]:

                if skeleton[yn, xn]:
                    G.add_edge((y, x), (yn, xn))

start = tuple(endpoint_coords[0])
end = tuple(endpoint_coords[1])

midline_path = nx.shortest_path(G, source = start, target = end)
print(f"This midline contains {len(midline_path)} units.")

# Midline spline
def fit_midline_spline(path, smoothing = 0, spacing = 4.0):
    path = np.array(path)
    x = path[:, 1]
    y = path[:, 0]

    dx = np.diff(x)
    dy = np.diff(y)

    distances = np.sqrt(dx**2 + dy**2)
    arc_length = np.concatenate([[0], np.cumsum(distances)])
    arc_length_norm = arc_length / arc_length[-1]

    tck, _ = splprep([x, y], u = arc_length_norm, s = smoothing)

    u_dense = np.linspace(0, 1, 1000)
    x_dense, y_dense = splev(u_dense, tck)
    points = np.column_stack((x_dense, y_dense))
    diffs = np.diff(points, axis = 0)
    step_lengths = np.linalg.norm(diffs, axis = 1)
    arc_lengths = np.concatenate([[0], np.cumsum(step_lengths)])
    total_length = arc_lengths[-1]
    
    num_samples = max(int(total_length / spacing), 2)
    target_lengths = np.linspace(0, total_length, num_samples)
    x_resampled = np.interp(target_lengths, arc_lengths, x_dense)
    y_resampled = np.interp(target_lengths, arc_lengths, y_dense)

    return list(zip(np.round(y_resampled).astype(int), np.round(x_resampled).astype(int)))

spline_path = fit_midline_spline(midline_path, smoothing = 5, spacing = 3)

# Extending spline ends to meet mask edge
def extend_to_mask_edge(mask, point, direction, max_steps = 10):

    y0, x0 = point
    dy, dx = direction
    h, w = mask.shape

    for step in range(1, max_steps +1):
        y = int(round(y0 + dy * step))
        x = int(round(x0 + dx * step))
        if not (0 <= y < h and 0 <= x < w):
            break
        if mask[y, x] == 0:
            y = int(round(y0 + dy * (step - 1)))
            x = int(round(x0 + dx * (step - 1)))
            return (y, x)
    return None

def interpolate_extension(p1, p2, spacing = 3):
    p1 = np.array(p1, dtype=float)
    p2 = np.array(p2, dtype=float)
    extension_length = np.linalg.norm(p2 - p1)
    num_points = max(int(np.round(extension_length / spacing)), 1)

    return [tuple(p1 + (p2 - p1) * t)
            for t in np.linspace(0, 1, num_points + 1)[1:]]

def deduplicate_rounded_path(path):
    seen = set()
    cleaned = []
    for pt in path:
        rounded = tuple(np.round(pt).astype(int))
        if rounded not in seen:
            seen.add(rounded)
            cleaned.append(pt)
    return cleaned

def get_unit_direction(p1, p2):
    dy = p2[0] - p1[0]
    dx = p2[1] - p1[1]
    norm = np.hypot(dy, dx)
    return (dy / norm, dx / norm) if norm != 0 else (0, 0)

def compute_normals(path, offset = 2):
    path = np.array(path).astype(float)
    N = len(path)
    tangents = np.zeros_like(path)

    for i in range(N):
        i1 = max(0, i - offset)
        i2 = min(N - 1, i + offset)
        delta = path[i2] - path[i1]
        norm = np.linalg.norm(delta)
        if norm > 1e-6:
            tangents[i] = delta / norm
        else:
            tangents[i] = [0, 0]

    normals = np.zeros_like(tangents)
    normals[:, 0] = -tangents[:, 1]         # dy = -dx
    normals[:, 1] = tangents[:, 0]          # dx = dy

    return normals

def compute_arc_lengths(path):
    path = np.array(path)
    diffs = np.diff(path, axis = 0)
    step_lengths = np.sqrt((diffs**2).sum(axis = 1))
    arc_lengths = np.concatenate([[0], np.cumsum(step_lengths)])

    return arc_lengths

def select_evenly_spaced_indices(arc_lengths, num_samples=10):
    
    arc_lengths = np.array(arc_lengths)
    target_positions = np.linspace(0, arc_lengths[-1], num_samples)
    indices = np.searchsorted(arc_lengths, target_positions)
    indices = sorted(set(indices))

    return indices

start_dir = get_unit_direction(spline_path[1], spline_path[0])
end_dir = get_unit_direction(spline_path[-2], spline_path[-1])

start_pad = extend_to_mask_edge(cleaned, spline_path[0], start_dir, max_steps = 20)
end_pad = extend_to_mask_edge(cleaned, spline_path[-1], end_dir, max_steps = 20)

spacing = 3.0

padded_spline = []

if start_pad:
    start_interp = interpolate_extension(spline_path[0], start_pad, spacing=spacing)
    padded_spline = start_interp + padded_spline

padded_spline += spline_path

if end_pad:
    end_interp = interpolate_extension(spline_path[-1], end_pad, spacing=spacing)
    padded_spline.extend(end_interp)

print("Corrected start padded point:", padded_spline[0])
print("Original spline_path[0]:", spline_path[0])

padded_spline = deduplicate_rounded_path(padded_spline)
print("First padded point (start):", padded_spline[0])
print("Second padded point:", padded_spline[1])
print("Start pad used:", start_pad)
print("First spline point:", spline_path[0])

padded_normals = compute_normals(padded_spline, offset=12)

arc_lengths = compute_arc_lengths(spline_path)
print(f"Total arc length: {arc_lengths[-1]:.2f} pixels")

if len(set(padded_spline)) != len(padded_spline):
    print("⚠️ Duplicate point detected in padded_spline")

padded_arcs = compute_arc_lengths(padded_spline)

# Width sampling
def sample_widths(mask, midline, normals, max_distance = 40, step_size = 0.5):
    mask = (mask > 0).astype(np.uint8)
    heights, widths_img = mask.shape

    widths = []
    edge_points = []

    for (y0, x0), (dy, dx) in zip(midline, normals):
        # Normalize direction
        norm = np.hypot(dy, dx)
        if norm == 0:
            widths.append(np.nan)
            edge_points.append((None, None))
            continue
        dy /= norm
        dx /= norm

        found1, found2 = False, False
        p1, p2 = None, None

        for d in np.arange(0, max_distance, step_size):
            y_pos = y0 + dy * d
            x_pos = x0 + dx * d
            y_neg = y0 - dy * d
            x_neg = x0 - dx * d

            # Bilinear sampling
            val_pos = map_coordinates(mask, [[y_pos], [x_pos]], order=1, mode='constant')[0]
            val_neg = map_coordinates(mask, [[y_neg], [x_neg]], order=1, mode='constant')[0]

            if not found1 and val_pos < 0.5:
                # Just crossed from inside to outside → use previous point
                y_back = y0 + dy * (d - step_size)
                x_back = x0 + dx * (d - step_size)
                p1 = (y_back, x_back)
                found1 = True

            if not found2 and val_neg < 0.5:
                y_back = y0 - dy * (d - step_size)
                x_back = x0 - dx * (d - step_size)
                p2 = (y_back, x_back)
                found2 = True

            if found1 and found2:
                break

        if found1 and found2:
            dist = np.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)
            widths.append(dist)
            edge_points.append((p1, p2))
        else:
            widths.append(np.nan)
            edge_points.append((None, None))

    return widths, edge_points        

widths, edge_pts = sample_widths(cleaned, padded_spline, padded_normals, max_distance = 50, step_size = 0.5)

# Volume

valid = [i for i, w in enumerate(widths) if not np.isnan(w)]
valid_widths = [widths[i] for i in valid]
valid_arcs = [padded_arcs[i] for i in valid]

def estimate_volume_ravg(widths, arc_position):
    volume = 0
    for i in range(len(widths) -1):
        w1, w2 = widths[1], widths[i+1]
        dx = arc_position[i + 1] - arc_position[i]
        r_avg = (w1 + w2 / 4)**2                        # divide by 4 bc we want r not 2r (= w)
        slice_volume = np.pi * r_avg * dx
        volume += slice_volume
    return volume

rotifer_volume_ravg = estimate_volume_ravg(valid_widths, valid_arcs)
print(f"Estimated rotifer volume using avg r: {rotifer_volume_ravg:.2f} pixels^3")

def estimate_volume_frustum(widths, arc_positions):
    volume = 0
    for i in range(len(widths) - 1):
        r1, r2 = widths[i] / 2, widths[i + 1] / 2
        dx = arc_positions[i + 1] - arc_positions[i]
        slice_vol = (np.pi * dx / 3) * (r1**2 + r1 * r2**2 + r2)
        volume += slice_vol
    return volume

rotifer_volume_frustum = estimate_volume_frustum(valid_widths, valid_arcs)
print(f"Estimated rotifer volume using frustum: {rotifer_volume_frustum:.2f} pixels^3")

print(f"{np.sum(np.isnan(widths))} out of {len(widths)} width samples were NaN")

# plot align check

print(len(padded_spline), len(padded_normals), len(widths), len(padded_arcs))

min_len = min(len(padded_spline), len(padded_normals), len(widths), len(padded_arcs))
padded_spline = padded_spline[:min_len]
padded_normals = padded_normals[:min_len]
widths = widths[:min_len]
padded_arcs = padded_arcs[:min_len]

print("First 5 arc_lengths:", padded_arcs[:5])
print("First 5 widths:", widths[:5])



# 3d wireframe comparison of vol calc method

def set_axes_equal(ax):
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = x_limits[1] - x_limits[0]
    y_range = y_limits[1] - y_limits[0]
    z_range = z_limits[1] - z_limits[0]

    max_range = max(x_range, y_range, z_range)

    x_middle = np.mean(x_limits)
    y_middle = np.mean(y_limits)
    z_middle = np.mean(z_limits)

    ax.set_xlim3d([x_middle - max_range/2, x_middle + max_range/2])
    ax.set_ylim3d([y_middle - max_range/2, y_middle + max_range/2])
    ax.set_zlim3d([z_middle - max_range/2, z_middle + max_range/2])

def build_local_frames(centers, initial_normal = np.array([1.0, 0.0, 0.0])):
    frames = []

    T0 = centers[1] - centers[0]
    T0 /= np.linalg.norm(T0)

    # Ensure N0 is orthogonal to T0
    T0_xy = T0[:2]
    perp_2d = np.array([-T0_xy[1], T0_xy[0]])
    perp_2d /= np.linalg.norm(perp_2d)

    # Embed this in 3D: keep Z = 0
    N0 = np.array([perp_2d[0], perp_2d[1], 0.0])
    B0 = np.cross(T0, N0)
    frames.append((T0, N0, B0))

    for i in range(1, len(centers)):
        T_prev = frames[-1][0]
        T = centers[i] - centers[i - 1]
        T /= np.linalg.norm(T)

        v = np.cross(T_prev, T)
        s = np.linalg.norm(v)
        c = np.dot(T_prev, T)

        if s < 1e-6:
            frames.append((T, *frames[-1][1:]))
            continue

        v /= s
        K = np.array([[0, -v[2], v[1]],
                      [v[2], 0, -v[0]],
                      [-v[1], v[0], 0]])

        R = np.eye(3) + K + K @ K * ((1 - c) / (s ** 2))

        N_prev, B_prev = frames[-1][1:]
        N = R @ N_prev
        N /= np.linalg.norm(N)
        B = R @ B_prev
        B /= np.linalg.norm(B)

        frames.append((T, N, B))

    return frames

def plot_oriented_3d_wireframe(spline_coords, arc_lengths, widths, step=2, n_points=40):

    fig = go.Figure()

    spline_coords = np.array(spline_coords).astype(float)
    arc_lengths = np.array(arc_lengths)
    centers = np.column_stack((spline_coords[:, 1], spline_coords[:, 0], arc_lengths))

    # Build initial tangent
    T0 = centers[1] - centers[0]
    T0 /= np.linalg.norm(T0)

    # Auto-orient first normal in image plane, orthogonal to T0
    T0_xy = T0[:2]
    perp_2d = np.array([-T0_xy[1], T0_xy[0]])
    perp_2d /= np.linalg.norm(perp_2d)
    N0 = np.array([perp_2d[0], perp_2d[1], 0.0])
    N0 -= T0 * np.dot(T0, N0)
    N0 /= np.linalg.norm(N0)
    B0 = np.cross(T0, N0)

    # Initialize frame list
    frames = [(T0, N0, B0)]

    # Propagate frame using parallel transport
    for i in range(1, len(centers)):
        T_prev = frames[-1][0]
        T = centers[i] - centers[i - 1]
        T /= np.linalg.norm(T)

        v = np.cross(T_prev, T)
        s = np.linalg.norm(v)
        c = np.dot(T_prev, T)

        if s < 1e-6:
            frames.append((T, *frames[-1][1:]))
            continue

        v /= s
        K = np.array([[0, -v[2], v[1]],
                      [v[2], 0, -v[0]],
                      [-v[1], v[0], 0]])
        R = np.eye(3) + K + K @ K * ((1 - c) / (s**2))

        N_prev, B_prev = frames[-1][1:]
        N = R @ N_prev
        B = R @ B_prev
        frames.append((T, N / np.linalg.norm(N), B / np.linalg.norm(B)))

    # Draw rings
    for i in range(0, len(widths), step):
        if np.isnan(widths[i]):
            continue

        center = centers[i]
        T, N, B = frames[i]
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
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Arc Length (Z)',
            aspectmode='data'
        ),
        title='3D Rotifer Reconstruction (Parallel Transport Frame)',
        margin=dict(l=0, r=0, t=40, b=0)
    )

    fig.show()


plot_oriented_3d_wireframe(padded_spline, padded_arcs, widths)


def plot_3d_rings(midline, widths, arcs, step=2, n_points=40):
    fig = go.Figure()

    for i in range(0, len(widths), step):
        if np.isnan(widths[i]):
            continue

        y, x = midline[i]
        z = arcs[i]
        r = widths[i] / 2

        theta = np.linspace(0, 2 * np.pi, n_points)
        circle_x = x + r * np.cos(theta)
        circle_y = y + r * np.sin(theta)
        circle_z = np.full_like(circle_x, z)

        fig.add_trace(go.Scatter3d(
            x = circle_x, y = circle_y, z = circle_z,
            mode = 'lines',
            line = dict(color = 'blue', width = 2),
            showlegend = False
        ))

    fig.update_layout(
        scene=dict(
            xaxis = dict(title = 'X'),
            yaxis = dict(title = 'Y'),
            zaxis = dict(title = 'Z'),
            aspectmode = 'data'
        ),
        title = '3D Wireframe',
        margin=dict(l=0, r=0, t=40, b=0)
    )

    fig.show()

# Composite overlay
overlay = np.zeros((*cleaned.shape, 3), dtype = np.uint8)
overlay[..., 0] = (cleaned == 0) * 255
overlay[..., 1] = cleaned * 255

for y0, x0 in midline_path:
    overlay[y0, x0, 2] = 255

for (y0, x0), (y1, x1) in zip(padded_spline[:-1], padded_spline[1:]):
    rr, cc = line(int(round(y0)), int(round(x0)), int(round(y1)), int(round(x1)))
    overlay[rr, cc, 0] = 0
    overlay[rr, cc, 1] = 0
    overlay[rr, cc, 2] = 255

for i, ((y, x), width) in enumerate(zip(padded_spline, widths)):
        if np.isnan(width):
            continue

# Plot
fig, axs = plt.subplots(2, 2, figsize = (10, 8))
axs = axs.flatten()
axs[0].imshow(cleaned, cmap = 'gray')
axs[0].set_title("Original Bionary ROI")
axs[0].axis('off')

axs[1].imshow(cleaned, cmap = 'gray')
axs[1].set_title("Cleaned Mask")
axs[1].axis('off')

axs[2].imshow(skeleton, cmap = 'grey')
axs[2].set_title("Skeleton")
axs[2].axis('off')

for i, (p1, p2) in enumerate(edge_pts):
    if p1 is not None and p2 is not None:
        rr, cc = line(int(round(p1[0])), int(round(p1[1])), int(round(p2[0])), int(round(p2[1])))
        overlay[rr, cc, :] = [255, 0, 0]


axs[3].imshow(overlay)
axs[3].set_title("Composite")
axs[3].axis('off')

for (y, x), (dy, dx) in zip(padded_spline[::1], padded_normals [::1]):     # plot all norms
    axs[3].arrow(x, y, dx*5, dy*5, color = 'cyan', head_width = 1)

plt.tight_layout()
plt.show()