#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <qasm_file> <device_json>" >&2
    echo "Env overrides:" >&2
    echo "  QASMTRANS_PY   (default: $REPO_ROOT/venv/bin/python)" >&2
    echo "  QISKIT_PY      (default: $REPO_ROOT/venv-qiskit313/bin/python)" >&2
    echo "  OPT_LEVEL      (default: 3)" >&2
    echo "  MODE           (default: ibmq)" >&2
    echo "  DISABLE_MAPOMATIC (default: 0)" >&2
    echo "  OPTIMIZE_1Q (default: 1)" >&2
    echo "  OPTIMIZE_2Q_CANCEL (default: 1)" >&2
    echo "  OPTIMIZE_COMMUTE_2Q (default: 1)" >&2
    echo "  OPTIMIZE_2Q_SYNTH (default: 1)" >&2
    echo "  ROUTING_DECAY (default: 0)" >&2
    echo "  ROUTING_TRIALS (default: 1)" >&2
    echo "  FAST_QUALITY (default: 1)" >&2
    echo "  FAST_QUALITY_EMBEDDINGS (default: 50)" >&2
    echo "  SEED (default: -1)" >&2
    exit 1
fi

QASM_FILE="$(readlink -f "$1")"
DEVICE_JSON="$(readlink -f "$2")"

QASMTRANS_PY="${QASMTRANS_PY:-$REPO_ROOT/venv/bin/python}"
QISKIT_PY="${QISKIT_PY:-$REPO_ROOT/venv-qiskit313/bin/python}"
OPT_LEVEL="${OPT_LEVEL:-3}"
MODE="${MODE:-ibmq}"
DISABLE_MAPOMATIC="${DISABLE_MAPOMATIC:-0}"
OPTIMIZE_1Q="${OPTIMIZE_1Q:-1}"
OPTIMIZE_2Q_CANCEL="${OPTIMIZE_2Q_CANCEL:-1}"
OPTIMIZE_COMMUTE_2Q="${OPTIMIZE_COMMUTE_2Q:-1}"
OPTIMIZE_2Q_SYNTH="${OPTIMIZE_2Q_SYNTH:-1}"
ROUTING_DECAY="${ROUTING_DECAY:-0}"
ROUTING_TRIALS="${ROUTING_TRIALS:-1}"
FAST_QUALITY="${FAST_QUALITY:-1}"
FAST_QUALITY_EMBEDDINGS="${FAST_QUALITY_EMBEDDINGS:-50}"
SEED="${SEED:-1}"

OUT_DIR="${OUT_DIR:-$REPO_ROOT/tmp/qasmtrans_compare}"
mkdir -p "$OUT_DIR"
QASMTRANS_PATH_FILE="${QASMTRANS_PATH_FILE:-$OUT_DIR/qasmtrans_path.txt}"

if [ ! -f "$QASM_FILE" ]; then
    echo "QASM file not found: $QASM_FILE" >&2
    exit 1
fi
if [ ! -f "$DEVICE_JSON" ]; then
    echo "Device JSON not found: $DEVICE_JSON" >&2
    exit 1
fi

export PYTHONPATH="$REPO_ROOT/python${PYTHONPATH:+:$PYTHONPATH}"

echo "QASM file: $QASM_FILE"
echo "Device JSON: $DEVICE_JSON"
echo

QASM_PATH="$QASM_FILE" \
DEVICE_JSON="$DEVICE_JSON" \
OUT_DIR="$OUT_DIR" \
QASMTRANS_PATH_FILE="$QASMTRANS_PATH_FILE" \
MODE="$MODE" \
DISABLE_MAPOMATIC="$DISABLE_MAPOMATIC" \
OPTIMIZE_1Q="$OPTIMIZE_1Q" \
OPTIMIZE_2Q_CANCEL="$OPTIMIZE_2Q_CANCEL" \
OPTIMIZE_COMMUTE_2Q="$OPTIMIZE_COMMUTE_2Q" \
OPTIMIZE_2Q_SYNTH="$OPTIMIZE_2Q_SYNTH" \
ROUTING_DECAY="$ROUTING_DECAY" \
ROUTING_TRIALS="$ROUTING_TRIALS" \
FAST_QUALITY="$FAST_QUALITY" \
FAST_QUALITY_EMBEDDINGS="$FAST_QUALITY_EMBEDDINGS" \
SEED="$SEED" \
"$QASMTRANS_PY" - <<'PY'
import os
import re
import time
from pathlib import Path

import re

import qasmtrans

qasm_path = Path(os.environ["QASM_PATH"])
device_json = Path(os.environ["DEVICE_JSON"])

opts = qasmtrans.TranspileOptions()
opts.mode = os.environ.get("MODE", "ibmq")
opts.backend_config = str(device_json)
opts.output_path = os.environ.get("OUT_DIR", "")
opts.disable_mapomatic = bool(int(os.environ.get("DISABLE_MAPOMATIC", "0")))
opts.optimize_1q = bool(int(os.environ.get("OPTIMIZE_1Q", "1")))
opts.optimize_2q_cancel = bool(int(os.environ.get("OPTIMIZE_2Q_CANCEL", "1")))
opts.optimize_commute_2q = bool(int(os.environ.get("OPTIMIZE_COMMUTE_2Q", "1")))
opts.optimize_2q_synth = bool(int(os.environ.get("OPTIMIZE_2Q_SYNTH", "1")))
opts.routing_decay = bool(int(os.environ.get("ROUTING_DECAY", "0")))
opts.routing_trials = int(os.environ.get("ROUTING_TRIALS", "1"))
opts.fast_quality_routing = bool(int(os.environ.get("FAST_QUALITY", "0")))
opts.fast_quality_max_embeddings = int(os.environ.get("FAST_QUALITY_EMBEDDINGS", "50"))
opts.seed = int(os.environ.get("SEED", "-1"))
opts.verbose = 1

qasm_text = qasm_path.read_text(encoding="utf-8")
start = time.perf_counter()
res = qasmtrans.transpile_qasm(qasm_text, opts)
wall_ms = (time.perf_counter() - start) * 1000.0

log = res.log if isinstance(res.log, str) else ""
time_match = re.search(r"total QASMTrans time:\\s*([-+]?\\d*\\.?\\d+)", log)
reported_ms = time_match.group(1) if time_match else "n/a"

out_qasm = res.output_qasm or ""
out_qasm_path = res.output_qasm_path or ""
if not out_qasm and res.output_qasm_path:
    try:
        out_qasm = Path(res.output_qasm_path).read_text(encoding="utf-8")
    except FileNotFoundError:
        out_qasm = ""
if not out_qasm:
    try:
        candidates = sorted(Path(opts.output_path).glob("transpiled_*.qasm"), key=lambda p: p.stat().st_mtime, reverse=True)
        if candidates:
            out_qasm = candidates[0].read_text(encoding="utf-8")
            out_qasm_path = str(candidates[0])
    except Exception:
        out_qasm = ""
if os.environ.get("QASMTRANS_PATH_FILE") and out_qasm_path:
    try:
        Path(os.environ["QASMTRANS_PATH_FILE"]).write_text(out_qasm_path, encoding="utf-8")
    except Exception:
        pass

if os.environ.get("DEBUG_QASMTRANS"):
    print(f"  debug_out_qasm_len={len(out_qasm)}")
    if out_qasm:
        first_line = out_qasm.splitlines()[0] if out_qasm.splitlines() else ""
        print(f"  debug_first_line={first_line}")
        for idx, line in enumerate(out_qasm.splitlines()[:5], start=1):
            print(f"  debug_line{idx}={line}")

print("QASMTrans:")
print(f"  wall_ms={wall_ms:.3f}")
print(f"  reported_ms={reported_ms}")
print(f"  output_qasm={out_qasm_path}")
PY

OPT_LEVEL="$OPT_LEVEL" \
QASM_PATH="$QASM_FILE" \
DEVICE_JSON="$DEVICE_JSON" \
QASMTRANS_PATH_FILE="$QASMTRANS_PATH_FILE" \
"$QISKIT_PY" - <<'PY'
import json
import os
import re
import time
from pathlib import Path

from qiskit import QuantumCircuit, transpile
from qiskit.transpiler import CouplingMap
try:
    from qiskit.qasm2 import dumps as qasm2_dumps
    from qiskit import qasm2
except Exception:  # pragma: no cover - fallback for older qiskit
    qasm2_dumps = None
    qasm2 = None

qasm_path = Path(os.environ["QASM_PATH"])
device_json = Path(os.environ["DEVICE_JSON"])
opt_level = int(os.environ.get("OPT_LEVEL", "3"))
qasmtrans_path_file = Path(os.environ.get("QASMTRANS_PATH_FILE", ""))

cfg = json.loads(device_json.read_text(encoding="utf-8"))
coupling = cfg.get("cx_coupling") or cfg.get("coupling_map")
if not coupling:
    raise SystemExit(f"No coupling_map/cx_coupling found in {device_json}")

edges = []
for entry in coupling:
    if isinstance(entry, str):
        parts = entry.replace("-", "_").split("_")
    else:
        parts = entry
    if len(parts) != 2:
        continue
    u, v = map(int, parts)
    edges.append((u, v))

basis = cfg.get("basis_gates") or ["rz", "sx", "x", "cx"]
basis_gates = list(dict.fromkeys(basis))

def ensure_ecr_definition(qasm_text: str) -> str:
    if "ecr" not in qasm_text.lower():
        return qasm_text
    if re.search(r"\\bgate\\s+ecr\\b", qasm_text, flags=re.IGNORECASE):
        return qasm_text
    gate_def = "gate ecr q0,q1 { s q0; sx q1; cx q0,q1; x q0; }"
    lines = qasm_text.splitlines()
    inserted = False
    out = []
    for line in lines:
        out.append(line)
        if not inserted and line.strip().lower().startswith("include"):
            out.append(gate_def)
            inserted = True
    if not inserted:
        out.insert(1 if lines else 0, gate_def)
    return "\\n".join(out) + "\\n"

def load_qasm_text(qasm_text: str) -> QuantumCircuit:
    qasm_text = ensure_ecr_definition(qasm_text)
    if qasm2 is not None:
        try:
            return qasm2.loads(qasm_text)
        except Exception:
            pass
    try:
        return QuantumCircuit.from_qasm_str(qasm_text)
    except Exception:
        tmp = Path("_tmp_qiskit_parse.qasm")
        tmp.write_text(qasm_text, encoding="utf-8")
        return QuantumCircuit.from_qasm_file(str(tmp))

def metrics_from_circuit(qc: QuantumCircuit):
    filtered = QuantumCircuit(qc.num_qubits, qc.num_clbits)
    for inst, qargs, cargs in qc.data:
        if inst.name in ("barrier", "measure", "delay", "reset"):
            continue
        filtered.append(inst, qargs, cargs)
    n1 = 0
    n2 = 0
    for inst, qargs, _ in filtered.data:
        if len(qargs) == 1:
            n1 += 1
        elif len(qargs) == 2:
            n2 += 1
    depth = filtered.depth()
    return n1, n2, depth

qc = QuantumCircuit.from_qasm_file(str(qasm_path))
start = time.perf_counter()
tc = transpile(
    qc,
    coupling_map=CouplingMap(edges),
    basis_gates=basis_gates,
    optimization_level=opt_level,
)
wall_ms = (time.perf_counter() - start) * 1000.0

qt_n1 = qt_n2 = qt_depth = 0
if qasmtrans_path_file and qasmtrans_path_file.exists():
    qt_path = qasmtrans_path_file.read_text(encoding="utf-8").strip()
    if qt_path:
        qt_text = Path(qt_path).read_text(encoding="utf-8")
        qt_circ = load_qasm_text(qt_text)
        qt_n1, qt_n2, qt_depth = metrics_from_circuit(qt_circ)

qk_n1, qk_n2, qk_depth = metrics_from_circuit(tc)

print("Qiskit:")
print(f"  opt_level={opt_level}")
print(f"  wall_ms={wall_ms:.3f}")
print(f"  depth={qk_depth} 1q={qk_n1} 2q={qk_n2}")
if qasmtrans_path_file and qasmtrans_path_file.exists():
    print("QASMTrans metrics (via Qiskit):")
    print(f"  depth={qt_depth} 1q={qt_n1} 2q={qt_n2}")
PY
