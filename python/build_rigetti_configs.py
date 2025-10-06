#!/usr/bin/env python3
"""Generate Rigetti Ankaa-3 device and pulse template configs from calibration data."""
from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

CALIBRATION_PATH = Path("data/devices/rigetti_calibrated_rx_iswap_iq.json")
DEVICE_OUTPUT = Path("data/devices/rigetti_ankaa3_device.json")
PULSE_OUTPUT = Path("data/devices/rigetti_ankaa3_pulses.json")

# Connectivity map provided by the calibration package (original qubit indices).
CONNECTIVITY: Dict[int, List[int]] = {
    0: [1, 7],
    1: [0, 2, 8],
    2: [1, 3, 9],
    3: [2, 4, 10],
    4: [3, 5, 11],
    5: [4, 6, 12],
    6: [5, 13],
    7: [0, 8, 14],
    8: [1, 7, 9, 15],
    9: [2, 8, 10, 16],
    10: [3, 9, 11],
    11: [4, 10, 12, 18],
    12: [5, 11, 13, 19],
    13: [6, 12, 20],
    14: [7, 15, 21],
    15: [8, 14, 16, 22],
    16: [9, 15, 17, 23],
    17: [16, 18, 24],
    18: [11, 17, 19, 25],
    19: [12, 18, 20, 26],
    20: [13, 19, 27],
    21: [14, 22, 28],
    22: [15, 21, 23, 29],
    23: [16, 22, 24, 30],
    24: [17, 23, 25, 31],
    25: [18, 24, 26, 32],
    26: [19, 25, 27, 33],
    27: [20, 26],
    28: [21, 29, 35],
    29: [22, 28, 30, 36],
    30: [23, 29, 31, 37],
    31: [24, 30, 38],
    32: [25, 33, 39],
    33: [26, 32, 34, 40],
    34: [33, 41],
    35: [28, 36],
    36: [29, 35, 37, 43],
    37: [30, 36, 38, 44],
    38: [31, 37, 39, 45],
    39: [32, 38, 40, 46],
    40: [33, 39, 41, 47],
    41: [34, 40],
    43: [36, 50],
    44: [37, 45, 51],
    45: [38, 44, 46, 52],
    46: [39, 45, 47, 53],
    47: [40, 46, 54],
    49: [50, 56],
    50: [43, 49, 51, 57],
    51: [44, 50, 52, 58],
    52: [45, 51, 53, 59],
    53: [46, 52, 54, 60],
    54: [47, 53, 55, 61],
    55: [54, 62],
    56: [49, 57, 63],
    57: [50, 56, 58, 64],
    58: [51, 57, 59, 65],
    59: [52, 58, 60, 66],
    60: [53, 59, 61, 67],
    61: [54, 60, 62, 68],
    62: [55, 61, 69],
    63: [56, 64, 70],
    64: [57, 63, 65, 71],
    65: [58, 64, 66, 72],
    66: [59, 65, 73],
    67: [60, 68, 74],
    68: [61, 67, 69, 75],
    69: [62, 68, 76],
    70: [63, 71, 77],
    71: [64, 70, 72, 78],
    72: [65, 71, 73, 79],
    73: [66, 72, 74, 80],
    74: [67, 73, 75, 81],
    75: [68, 74, 76, 82],
    76: [69, 75, 83],
    77: [70, 78],
    78: [71, 77, 79],
    79: [72, 78, 80],
    80: [73, 79, 81],
    81: [74, 80, 82],
    82: [75, 81, 83],
    83: [76, 82],
}

PI_OVER_TWO = 1.5707963267948966
ANGLE_TOL = 1e-9


def angle_matches(value: float, target: float, tol: float = ANGLE_TOL) -> bool:
    return math.isclose(value, target, rel_tol=0.0, abs_tol=tol)


def angle_label(angle: float) -> str:
    if angle_matches(angle, PI_OVER_TWO):
        return "pi_over_2"
    if angle_matches(angle, -PI_OVER_TWO):
        return "neg_pi_over_2"
    if angle_matches(angle, math.pi):
        return "pi"
    if angle_matches(angle, -math.pi):
        return "neg_pi"
    formatted = f"{angle:+.6f}".replace("+", "pos_").replace("-", "neg_")
    return formatted.replace(".", "_")


def angle_display(angle: float) -> str:
    if angle_matches(angle, PI_OVER_TWO):
        return "pi/2"
    if angle_matches(angle, -PI_OVER_TWO):
        return "-pi/2"
    if angle_matches(angle, math.pi):
        return "pi"
    if angle_matches(angle, -math.pi):
        return "-pi"
    return f"{angle:.6f}"


def load_calibration() -> dict:
    with CALIBRATION_PATH.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def build_qubit_mapping(connectivity: Dict[int, Iterable[int]]) -> Tuple[Dict[int, int], List[int]]:
    physical = set(connectivity.keys())
    for nbrs in connectivity.values():
        physical.update(nbrs)
    ordered = sorted(physical)
    mapping = {physical_id: idx for idx, physical_id in enumerate(ordered)}
    return mapping, ordered


def extract_rx_calibrations(calibrations: List[dict]) -> Dict[int, List[dict]]:
    rx_entries: Dict[int, List[dict]] = defaultdict(list)
    for entry in calibrations:
        if entry.get("gate") != "Rx":
            continue
        qubit = entry.get("qubits", [None])[0]
        if qubit is None:
            continue
        rx_entries[qubit].append(entry)
    return rx_entries


def extract_iswap_calibrations(calibrations: List[dict]) -> Dict[Tuple[int, int], dict]:
    iswap_entries: Dict[Tuple[int, int], dict] = {}
    for entry in calibrations:
        if entry.get("gate") != "ISwap":
            continue
        qubit_pair = entry.get("qubits", [])
        if len(qubit_pair) != 2:
            continue
        key = tuple(qubit_pair)
        iswap_entries[key] = entry
    return iswap_entries


def compute_duration(samples: List[float], dt: float) -> float:
    return float(len(samples)) * float(dt)


def compute_amplitude(samples_i: List[float], samples_q: List[float]) -> float:
    peak_i = max((abs(val) for val in samples_i), default=0.0)
    peak_q = max((abs(val) for val in samples_q), default=0.0)
    return float(max(peak_i, peak_q))


def generate_configs() -> None:
    calibration_doc = load_calibration()
    mapping, ordered_phys = build_qubit_mapping(CONNECTIVITY)

    rx_entries = extract_rx_calibrations(calibration_doc["calibrations"])
    iswap_entries = extract_iswap_calibrations(calibration_doc["calibrations"])

    missing_rx = set(mapping.keys()) - set(rx_entries.keys())
    if missing_rx:
        raise RuntimeError(f"Missing Rx calibrations for qubits: {sorted(missing_rx)}")

    pulse_definitions: List[dict] = []
    gate_lens: Dict[str, float] = {}
    gate_errs: Dict[str, float] = {}

    for physical_qubit, entries in rx_entries.items():
        mapped_qubit = mapping[physical_qubit]
        longest_rx_duration = 0.0
        for entry in entries:
            iq_data = entry.get("iq", {})
            samples_i = [float(val) for val in iq_data.get("i", [])]
            samples_q = [float(val) for val in iq_data.get("q", [])]
            dt = float(entry.get("dt", 1e-9))
            duration = compute_duration(samples_i, dt)
            amplitude = compute_amplitude(samples_i, samples_q)
            longest_rx_duration = max(longest_rx_duration, duration)
            angle = float(entry.get("params", {}).get("angle", 0.0))
            label = angle_label(angle)
            pulse_definitions.append(
                {
                    "id": f"rx_q{mapped_qubit}_{label}",
                    "gate": "rx",
                    "qubits": [mapped_qubit],
                    "shape": "arbitrary",
                    "waveform_type": "arbitrary",
                    "width": duration,
                    "amplitude": amplitude,
                    "samples_i": samples_i,
                    "samples_q": samples_q,
                    "parameters": {"theta": angle},
                    "note": (
                        f"source_qubit={physical_qubit}, theta={angle_display(angle)}"
                    ),
                }
            )
        pulse_definitions.append(
            {
                "id": f"rz_q{mapped_qubit}",
                "gate": "rz",
                "qubits": [mapped_qubit],
                "shape": "virtual",
                "waveform_type": "virtual",
                "width": 0.0,
                "amplitude": 0.0,
                "virtual": True,
                "note": f"virtual frame change for source_qubit={physical_qubit}",
            }
        )
        gate_lens[f"rx{mapped_qubit}"] = longest_rx_duration
        gate_errs[f"rx{mapped_qubit}"] = 0.0
        gate_lens[f"rz{mapped_qubit}"] = 0.0
        gate_errs[f"rz{mapped_qubit}"] = 0.0
        gate_lens[f"reset{mapped_qubit}"] = 0.0
        gate_errs[f"reset{mapped_qubit}"] = 0.0

    for key, entry in iswap_entries.items():
        physical_q0, physical_q1 = key
        mapped_q0 = mapping[physical_q0]
        mapped_q1 = mapping[physical_q1]
        iq_data = entry.get("iq", {})
        samples_i = [float(val) for val in iq_data.get("i", [])]
        samples_q = [float(val) for val in iq_data.get("q", [])]
        dt = float(entry.get("dt", 1e-9))
        duration = compute_duration(samples_i, dt)
        amplitude = compute_amplitude(samples_i, samples_q)
        pulse_definitions.append(
            {
                "id": f"iswap_q{mapped_q0}_q{mapped_q1}",
                "gate": "iswap",
                "qubits": [mapped_q0, mapped_q1],
                "shape": "arbitrary",
                "waveform_type": "arbitrary",
                "width": duration,
                "amplitude": amplitude,
                "samples_i": samples_i,
                "samples_q": samples_q,
                "note": f"source_pair=({physical_q0},{physical_q1})",
            }
        )
        gate_lens[f"iswap{mapped_q0}_{mapped_q1}"] = duration
        gate_errs[f"iswap{mapped_q0}_{mapped_q1}"] = 0.0

    cx_coupling = sorted(
        {
            f"{mapping[src]}_{mapping[dst]}"
            for src, nbrs in CONNECTIVITY.items()
            for dst in nbrs
        }
    )

    device_document = {
        "name": "rigetti_ankaa3_calibrated",
        "version": calibration_doc.get("retrieved_at_utc", "0.1"),
        "num_qubits": len(ordered_phys),
        "basis_gates": ["rx", "rz", "iswap", "reset"],
        "gate_lens": gate_lens,
        "gate_errs": gate_errs,
        "cx_coupling": cx_coupling,
        "physical_qubits": ordered_phys,
    }

    pulse_document = {
        "name": "rigetti_ankaa3_calibrated",
        "version": calibration_doc.get("retrieved_at_utc", "0.1"),
        "num_qubits": len(ordered_phys),
        "basis_gates": ["rx", "rz", "iswap"],
        "pulse_definitions": pulse_definitions,
    }

    DEVICE_OUTPUT.write_text(json.dumps(device_document, indent=2) + "\n", encoding="utf-8")
    PULSE_OUTPUT.write_text(json.dumps(pulse_document, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    generate_configs()
