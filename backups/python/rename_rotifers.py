import os
import re
from pathlib import Path
from PIL import Image

# === CONFIGURATION ===
input_root = Path(r"G:\Research\imageAnalysis\batch 1")  # <-- CHANGE THIS
output_root = Path(r"G:\Research\imageAnalysis\raw images")

# === PROCESSING ===
pattern = re.compile(r'^(.*?)\.MP4\..*?\.Still\d+[-_](\w+)\.png$', re.IGNORECASE)

for batch_folder in input_root.iterdir():
    if not batch_folder.is_dir():
        continue

    print(f"Processing batch: {batch_folder.name}")
    images = sorted(batch_folder.glob("*.png"))
    counter = 1

    for img_path in images:
        match = pattern.match(img_path.name)
        if not match:
            print(f"Skipping (no match): {img_path.name}")
            continue

        batch_id = match.group(1)
        posture = match.group(2)
        new_still_num = f"{counter:04d}"
        new_name = f"{batch_id}_{new_still_num}_{posture}.tif"

        # Create subfolder if it doesn't exist
        output_subfolder = output_root / batch_folder.name
        output_subfolder.mkdir(parents=True, exist_ok=True)

        new_path = output_subfolder / new_name

        # Open image and convert to TIFF
        with Image.open(img_path) as im:
            im.save(new_path, format="TIFF")


        counter += 1

    print(f"Finished batch '{batch_folder.name}' with {counter - 1} files.\n")
