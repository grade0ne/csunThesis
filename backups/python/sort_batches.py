import shutil
from pathlib import Path

# === CONFIGURATION ===
input_root = Path(r"G:\Research\imageAnalysis\resized images")  # <-- change this to your resized folder
segmentable_classes = {"s", "sf", "sm", "sfm", "w"}

# === SCRIPT ===
for batch_folder in input_root.iterdir():
    if not batch_folder.is_dir():
        continue

    segmentable_dir = batch_folder / "segmentable"
    manual_dir = batch_folder / "manual"
    segmentable_dir.mkdir(exist_ok=True)
    manual_dir.mkdir(exist_ok=True)

    for tif_file in batch_folder.glob("*.tif"):
        parts = tif_file.stem.split("_")
        if len(parts) < 3:
            print(f"Skipping (unreadable name): {tif_file.name}")
            continue

        posture = parts[-1].lower()
        dest_folder = segmentable_dir if posture in segmentable_classes else manual_dir
        dest_path = dest_folder / tif_file.name

        print(f"Moving {tif_file.name} → {dest_folder.name}")
        shutil.move(str(tif_file), str(dest_path))

print("✅ Done.")
