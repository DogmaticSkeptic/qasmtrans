#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="$REPO_ROOT/build"
OUTPUT_CSV="$REPO_ROOT/mapomatic_fidelity_comparison.csv"
OUTPUT_DIR="$REPO_ROOT/data/output_qasm_file"
DEVICE_JSON="$REPO_ROOT/data/devices/ibm_brisbane.json"
NWQSIM_EXE="${NWQSIM_EXE:-$REPO_ROOT/../NWQ-Sim/build/qasm/nwq_qasm}"
NWQSIM_BACKEND="${NWQSIM_BACKEND:-CPU}"
NWQSIM_SHOTS="${NWQSIM_SHOTS:-1024}"
NWQSIM_EXTRA_ARGS_STR="${NWQSIM_EXTRA_ARGS:-}"
declare -a NWQSIM_EXTRA_ARGS=()
if [[ -n "$NWQSIM_EXTRA_ARGS_STR" ]]; then
    read -r -a NWQSIM_EXTRA_ARGS <<<"$NWQSIM_EXTRA_ARGS_STR"
fi

declare -a inputs=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i)
            if [[ $# -lt 2 ]]; then
                echo "-i requires a path" >&2
                exit 1
            fi
            abs_path=$(python3 - "$2" <<'PY'
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
)
            inputs+=("$abs_path")
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

if [[ ${#inputs[@]} -eq 0 ]]; then
    echo "No input circuits provided (-i)" >&2
    exit 1
fi

if [[ ! -x $NWQSIM_EXE ]]; then
    echo "NWQ-Sim executable not found at $NWQSIM_EXE" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
cd "$BUILD_DIR"

echo "Circuit,Mode,PhysicalQubits,AverageFidelity" > "$OUTPUT_CSV"

measure_fidelity() {
    local device_json="$1"
    local circuit_qasm="$2"
    local cmd=("$NWQSIM_EXE" --backend "$NWQSIM_BACKEND" --shots "$NWQSIM_SHOTS" --sim dm --device "$device_json" -q "$circuit_qasm" --fidelity)
    if ((${#NWQSIM_EXTRA_ARGS[@]})); then
        cmd+=("${NWQSIM_EXTRA_ARGS[@]}")
    fi
    local log
    if ! log="$("${cmd[@]}" 2>&1)"; then
        printf '%s\n' "$log" >&2
        echo "NWQ-Sim execution failed" >&2
        exit 1
    fi
    local fidelity
    fidelity=$(python3 -c '
import re, sys
text = sys.stdin.read()
match = re.search(r"State Fidelity:\s*([0-9.eE+-]+)", text)
if not match:
    sys.exit(1)
print(match.group(1))
' <<<"$log") || {
        printf '%s\n' "$log" >&2
        echo "Failed to extract fidelity" >&2
        exit 1
    }
    printf '%s' "$fidelity"
}

run_case() {
    local circuit="$1"
    local circuit_name
    circuit_name="$(basename "$circuit" .qasm)"
    local output_base="$OUTPUT_DIR/${circuit_name}_${2}"
    local args=("./qasmtrans" "-i" "$circuit" "-m" "ibmq" "-c" "$DEVICE_JSON")
    args+=("-o" "${output_base}.qasm")
    if [[ "$2" == "off" ]]; then
        args+=("--disable_mapomatic")
    fi
    if ! "${args[@]}" > /dev/null; then
        echo "qasmtrans failed for $circuit ($2)" >&2
        exit 1
    fi
    local subchip_dir="${output_base}_subchips"
    local device_json="$subchip_dir/circuit00_subchip.json"
    local circuit_qasm="$subchip_dir/circuit00_subchip.qasm"
    local phys
    phys=$(python3 - "$device_json" <<'PY'
import json,sys
with open(sys.argv[1]) as f:
    data=json.load(f)
print(len(data.get('local_to_global', [])))
PY
) || phys=""
    local fidelity
    fidelity=$(measure_fidelity "$device_json" "$circuit_qasm")
    printf '%s,%s,%s,%s\n' "$(basename "$circuit")" "$2" "$phys" "$fidelity" >> "$OUTPUT_CSV"
}

for circuit in "${inputs[@]}"; do
    run_case "$circuit" "off"
    run_case "$circuit" "on"
done

echo "CSV written to $OUTPUT_CSV"
