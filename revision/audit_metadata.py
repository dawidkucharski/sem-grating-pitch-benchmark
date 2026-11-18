#!/usr/bin/env python3
"""Audit the JEOL sidecars without treating directory labels as magnifications."""

import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "revision" / "analysis"
GROUPS = (400, 500, 750, 1000)
DETECTORS = ("SED", "BED-S")


def read_sidecar(path):
    raw = path.read_bytes()
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = raw.decode("cp1252")
    values = {}
    for line in text.splitlines():
        if line.startswith("$"):
            key, _, value = line.partition(" ")
            values[key] = value.strip()
    return values, hashlib.sha256(raw).hexdigest()


def length_um(value):
    match = re.fullmatch(r"([0-9.eE+\-]+)\s*(um|\u00b5m|\u03bcm|nm|mm)", value)
    if not match:
        raise ValueError(f"Unsupported metadata length: {value!r}")
    number, unit = match.groups()
    factor = {"um": 1, "\u00b5m": 1, "\u03bcm": 1, "nm": 0.001, "mm": 1000}[unit]
    return float(number) * factor


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def audit():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    records = []
    for group in GROUPS:
        for detector in DETECTORS:
            folder = ROOT / f"mag{group}_SED-BED" / detector
            images = sorted(folder.glob("*.png"))
            sidecars = sorted(folder.glob("*.txt"))
            if {path.stem for path in images} != {path.stem for path in sidecars}:
                raise ValueError(f"Image/sidecar mismatch in {folder}")
            for path in sidecars:
                metadata, checksum = read_sidecar(path)
                width, height = map(int, metadata["$CM_IMAGE_SIZE"].split())
                fov_x, fov_y = map(length_um, metadata["$CM_FIELD_OF_VIEW"].split())
                stage = list(map(float, metadata["$CM_STAGE_POSITION"].split()))
                with Image.open(path.with_suffix(".png")) as image:
                    if image.size != (width, height):
                        raise ValueError(f"Image/metadata dimensions disagree: {path}")
                records.append({
                    "group_label": group,
                    "detector": detector,
                    "image_path": path.with_suffix(".png").relative_to(ROOT).as_posix(),
                    "sidecar_sha256": checksum,
                    "acquisition_magnification": float(metadata["$CM_MAG"]),
                    "width_px": width,
                    "height_px": height,
                    "fov_x_um": fov_x,
                    "fov_y_um": fov_y,
                    "scale_x_um_per_px": fov_x / width,
                    "scale_y_um_per_px": fov_y / height,
                    "stage_x_um": stage[0] * 1000,
                    "stage_y_um": stage[1] * 1000,
                    "stage_position_raw": metadata["$CM_STAGE_POSITION"],
                    "timestamp": metadata["$CM_DATE"] + "T" + metadata["$CM_TIME"],
                    "column_voltage_kV": float(metadata["$CM_ACCEL_VOLT"]),
                    "stage_bias_V": float(metadata["$SM_STAGE_BIAS_VOLT"]),
                    "working_distance_mm": float(metadata["$SM_WD"]),
                    "scan_time_setting": metadata["$SM_SCAN_TIME"],
                    "integration_number": metadata["$SM_INTEGRATION_NUMBER"],
                    "scan_rotation": metadata["$SM_SCAN_ROTATION"],
                })
    grouped = defaultdict(list)
    for record in records:
        grouped[(record["group_label"], record["detector"])].append(record)
    summaries = []
    for (group, detector), subset in sorted(grouped.items()):
        summary = {"group_label": group, "detector": detector, "images": len(subset)}
        for field in ("acquisition_magnification", "width_px", "height_px", "fov_x_um", "fov_y_um",
                      "scale_x_um_per_px", "scale_y_um_per_px", "column_voltage_kV", "stage_bias_V",
                      "scan_time_setting", "integration_number", "scan_rotation"):
            values = sorted({record[field] for record in subset})
            summary[field] = ";".join(map(str, values))
        for axis in ("x", "y"):
            values = [record[f"stage_{axis}_um"] for record in subset]
            summary[f"stage_{axis}_range_um"] = max(values) - min(values)
        summary["first_timestamp"] = min(record["timestamp"] for record in subset)
        summary["last_timestamp"] = max(record["timestamp"] for record in subset)
        summaries.append(summary)
    write_csv(OUTPUT / "sem_metadata.csv", records)
    write_csv(OUTPUT / "sem_metadata_summary.csv", summaries)
    print(json.dumps({
        "images_audited": len(records),
        "acquisition_magnifications": sorted({record["acquisition_magnification"] for record in records}),
        "horizontal_scales_um_per_px": sorted({record["scale_x_um_per_px"] for record in records}),
        "groups": summaries,
    }, indent=2))


if __name__ == "__main__":
    assert abs(length_um("17.07um") / 1280 - 0.0133359375) < 1e-14
    assert length_um("903000nm") == 903
    audit()