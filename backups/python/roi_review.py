import os
import cv2
import pandas as pd
import numpy as np
from roifile import ImagejRoi
from zipfile import ZipFile
from datetime import datetime

BASE_IMAGE_DIR = r"G:\Research\imageAnalysis\resized images"
BASE_ROI_DIR = r"G:\Research\imageAnalysis\partical analysis results\ROIs"
BASE_RESULT_DIR = r"G:\Research\imageAnalysis\partical analysis results"
OUTPUT_DIR = r"G:\Research\imageAnalysis\accepted results"
LOG_PATH = os.path.join(OUTPUT_DIR, f"review_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv")
os.makedirs(OUTPUT_DIR, exist_ok=True)

POSTURE_LABELS = {
    "s": "Stretched", "sf": "Stretched Feeding", "sm": "Stretched Multiple",
    "w": "Swimming", "wm": "Swimming Multiple"
}

KEY_HELP = "[Y]: Accept  [N]: Reject  [<]: Back  [>]: Fwd  [Tab] Toggle View"

def load_rois(zip_path):
    rois = []
    with ZipFile(zip_path, 'r') as archive:
        for name in archive.namelist():
            if name.lower().endswith('.roi'):
                with archive.open(name) as f:
                    data = f.read()
                    if data:
                        roi = ImagejRoi.frombytes(data)
                        rois.append((name, roi))
    return rois

def draw_roi_overlay(img, roi, label=None, status="undecided", show_overlay=True,
                     roi_index=0, total_rois=0, frame_id="", sample_id=""):
    overlay = img.copy()

    if show_overlay:
        if roi is not None:
            points = roi.coordinates()
            color = {"accepted": (0,255,0), "rejected": (0,0,255)}.get(status, (0,255,255))
            cv2.polylines(overlay, [points], isClosed=True, color=color, thickness=2)

        # Top-left info
        info_lines = [
            f"{sample_id}",
            f"{frame_id}",
            f"ROI: {roi_index+1}/{total_rois}"
        ]
        for i, line in enumerate(info_lines):
            y = 30 + i * 30
            cv2.putText(overlay, line, (15, y), cv2.FONT_HERSHEY_SIMPLEX,
                        0.8, (255,255,255), 1, lineType=cv2.LINE_AA)

        # Top-right posture label
        if label:
            font_scale = 1
            thickness = 2
            text_size = cv2.getTextSize(label, cv2.FONT_HERSHEY_SIMPLEX, font_scale, thickness)[0]
            x = overlay.shape[1] - text_size[0] - 15
            y = 40
            cv2.putText(overlay, label, (x, y), cv2.FONT_HERSHEY_SIMPLEX,
                        font_scale, (255,255,255), thickness, lineType=cv2.LINE_AA)

        # Bottom key help
        cv2.putText(overlay, KEY_HELP, (15, overlay.shape[0] - 15),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255,255,255), 1, lineType=cv2.LINE_AA)

    return overlay

def review_single_image(image_path, roi_zip_path, result_csv_path, output_path,
                        log_entries, previous_state=None):
    img = cv2.imread(image_path)
    rois = load_rois(roi_zip_path)
    results = pd.read_csv(result_csv_path)

    parts = os.path.splitext(os.path.basename(image_path))[0].split("_")
    sample_id = parts[0]
    frame_id = parts[1]
    posture_code = parts[2]
    posture_label = POSTURE_LABELS.get(posture_code.lower(), posture_code)

    decisions = previous_state if previous_state and len(previous_state) == len(rois) else [None] * len(rois)
    idx = 0
    show_overlay = True

    while 0 <= idx < len(rois):
        name, roi = rois[idx]
        status = decisions[idx] or "undecided"
        overlay = draw_roi_overlay(img, roi, posture_label, status, show_overlay,
                                   roi_index=idx, total_rois=len(rois),
                                   frame_id=frame_id, sample_id=sample_id)

        while True:
            cv2.imshow("ROI Review", overlay)
            key = cv2.waitKey(0)
            if key == ord('y'):
                decisions[idx] = "accepted"
                print(f"[✓] ROI {idx+1}/{len(rois)}: accepted")
            elif key == ord('n'):
                decisions[idx] = "rejected"
                print(f"[X] ROI {idx+1}/{len(rois)}: rejected")
            elif key == 9:  # Tab
                show_overlay = not show_overlay
            elif key == ord('.'):
                idx += 1
                break
            elif key == ord(','):
                if idx == 0:
                    cv2.destroyAllWindows()
                    return decisions, "back"
                else:
                    idx -= 1
                    break

            overlay = draw_roi_overlay(img, roi, posture_label, decisions[idx], show_overlay,
                                       roi_index=idx, total_rois=len(rois),
                                       frame_id=frame_id, sample_id=sample_id)

    cv2.destroyAllWindows()
    if any(d == "accepted" for d in decisions):
        results.iloc[[i for i, d in enumerate(decisions) if d == "accepted"]].to_csv(output_path, index=False)

    for i, d in enumerate(decisions):
        log_entries.append({
            "image": os.path.basename(image_path),
            "roi_index": i,
            "decision": d if d else "skipped",
            "posture": posture_label
        })

    return decisions, "next"

def batch_review():
    files = []
    for fname in sorted(os.listdir(BASE_ROI_DIR)):
        if fname.endswith("_rois.zip"):
            base = fname.replace("_rois.zip", "")
            parts = base.split("_")
            if len(parts) >= 3:
                sid, frame_id, posture = parts[0], parts[1], parts[2]
                img_path = os.path.join(BASE_IMAGE_DIR, sid, "segmentable", f"{sid}_{frame_id}_{posture}.tif")
                res_path = os.path.join(BASE_RESULT_DIR, base + "_results.csv")
                roi_path = os.path.join(BASE_ROI_DIR, fname)
                out_path = os.path.join(OUTPUT_DIR, base + "_accepted.csv")
                if os.path.exists(img_path) and os.path.exists(res_path):
                    files.append((img_path, roi_path, res_path, out_path))

    idx = 0
    memory = [None] * len(files)
    log_entries = []

    while 0 <= idx < len(files):
        img, roi, res, out = files[idx]
        prev_state = memory[idx]
        decisions, nav = review_single_image(img, roi, res, out, log_entries, prev_state)
        memory[idx] = decisions

        if nav == "back":
            if idx > 0:
                idx -= 1
            else:
                print("[!] Already at first image")
        else:
            idx += 1

    pd.DataFrame(log_entries).to_csv(LOG_PATH, index=False)
    print(f"[✓] Log saved to: {LOG_PATH}")

if __name__ == "__main__":
    batch_review()
