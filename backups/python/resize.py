import os
from pathlib import Path
from PIL import Image

# === CONFIGURATION ===
input_root = Path(r"G:\Research\imageAnalysis\raw images")  # <-- Your renamed image folder
output_root = Path(r"G:\Research\imageAnalysis\resized images")  # <-- Where to save resized versions

scale_factor = 0.5

# === PROCESSING ===
for batch_folder in input_root.iterdir():
    if not batch_folder.is_dir():
        continue

    output_batch = output_root / batch_folder.name
    output_batch.mkdir(parents=True, exist_ok=True)

    for tif_file in batch_folder.glob("*.tif"):
        with Image.open(tif_file) as im:
            new_size = (int(im.width * scale_factor), int(im.height * scale_factor))
            resized = im.resize(new_size, Image.BILINEAR)

            out_path = output_batch / tif_file.name
            resized.save(out_path, format="TIFF")

    print(f"Finished resizing: {batch_folder.name}")
