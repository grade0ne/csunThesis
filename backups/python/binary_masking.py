from pathlib import Path
import numpy as np
import tifffile as tiff

# === CONFIG ===
input_dir = Path(r"G:/Research/imageAnalysis/segmented images")
output_dir = Path(r"G:/Research/imageAnalysis/binary masks")
output_dir.mkdir(parents=True, exist_ok=True)

# === PROCESS ALL .tiff FILES ===
tif_files = list(input_dir.rglob("*.tiff"))
print(f"Found {len(tif_files)} files.")

for tif_path in tif_files:
    seg = tiff.imread(tif_path)
    
    # Create binary mask: rotifer class = 1.0
    rotifer_mask = (np.isclose(seg, 1.0)).astype(np.uint8) * 255
    
    out_path = output_dir / (tif_path.stem + "_BIN.tiff")
    tiff.imwrite(out_path, rotifer_mask, dtype=np.uint8)
    print(f"✅ Saved: {out_path.name}")