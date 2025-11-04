#!/usr/bin/env bash
# Collect fidelity and latency metrics for the 4-qubit benchmark algorithms.
# For each algorithm in data/test_benchmark/experiments/{qft4,shor4,adder4,w4,bv4,ghz4},
# the script performs an end-to-end workflow:
#   1. Compile the logical circuit with the baseline Rigetti templates to obtain a baseline schedule.
#   2. Optimise merged pulses with merge_pulse.py, producing per-algorithm merged templates/devices.
#   3. Recompile the circuit against the merged artefacts to obtain a post-merge schedule.
#   4. Simulate both schedules with simulate_ankaa9q.py.
#   5. Emit fidelity/latency metrics to data/output/alg4/alg4_metrics.csv.
#
# Requirements:
#   - jq
#   - QASMTrans build artefacts in ./build
#   - Conda environment "QASMtrans" available via /opt/homebrew/anaconda3/bin/activate

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"

if ! command -v jq >/dev/null 2>&1; then
    echo "Error: jq is required but not found in PATH." >&2
    exit 1
fi

CONDA_ACTIVATE="/opt/homebrew/anaconda3/bin/activate"
ENV_NAME="QASMtrans"
if [ ! -f "$CONDA_ACTIVATE" ]; then
    echo "Error: conda activate script not found at $CONDA_ACTIVATE" >&2
    exit 1
fi

# shellcheck disable=SC1091
source "$CONDA_ACTIVATE" "$ENV_NAME"

export MPLCONFIGDIR="${REPO_DIR}/.mplconfig"
mkdir -p "$MPLCONFIGDIR"

usage() {
    cat <<EOF
Usage: ${0##*/} [options]

Options:
  --merge-limit N        Maximum number of merged candidates to keep (default: env MERGE_LIMIT or 1)
  --base-device PATH     Baseline device JSON (default: env BASE_DEVICE or ${BASE_DEVICE_DEFAULT})
  --base-pulses PATH     Baseline pulse template JSON (default: env BASE_PULSE_TEMPLATE or ${BASE_PULSES_DEFAULT})
  --plot / --no-plot     Enable or disable PDF plot generation (default: enabled)
  -h, --help             Show this help message and exit
EOF
}

get_qasm_stub() {
    case "$1" in
        qft)  printf 'qft4' ;;
        shor) printf 'shor4' ;;
        adder) printf 'adder4' ;;
        w)    printf 'w4' ;;
        bv)   printf 'bv4' ;;
        ghz)  printf 'ghz4' ;;
        *)    return 1 ;;
    esac
}

resolve_path() {
    local candidate=$1
    if [[ $candidate = /* ]]; then
        printf '%s\n' "$candidate"
    else
        printf '%s\n' "${REPO_DIR}/${candidate}"
    fi
}

BASE_DEVICE_DEFAULT="data/devices/rigetti_ankaa4q_ideal_device.json"
BASE_PULSES_DEFAULT="data/devices/rigetti_ankaa4q_ideal_pulses.json"

BASE_DEVICE_INPUT="${BASE_DEVICE:-$BASE_DEVICE_DEFAULT}"
BASE_PULSES_INPUT="${BASE_PULSE_TEMPLATE:-$BASE_PULSES_DEFAULT}"
MERGE_LIMIT_INPUT="${MERGE_LIMIT:-1}"
DO_PLOT=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --merge-limit)
            [[ $# -ge 2 ]] || { echo "Error: --merge-limit requires a value." >&2; usage; exit 1; }
            MERGE_LIMIT_INPUT="$2"
            shift 2
            ;;
        --base-device)
            [[ $# -ge 2 ]] || { echo "Error: --base-device requires a value." >&2; usage; exit 1; }
            BASE_DEVICE_INPUT="$2"
            shift 2
            ;;
        --base-pulses)
            [[ $# -ge 2 ]] || { echo "Error: --base-pulses requires a value." >&2; usage; exit 1; }
            BASE_PULSES_INPUT="$2"
            shift 2
            ;;
        --plot)
            DO_PLOT=1
            shift
            ;;
        --no-plot)
            DO_PLOT=0
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: unknown option '$1'." >&2
            usage
            exit 1
            ;;
    esac
done

BASE_DEVICE_PATH="$(resolve_path "$BASE_DEVICE_INPUT")"
BASE_PULSES_PATH="$(resolve_path "$BASE_PULSES_INPUT")"

if [ ! -f "$BASE_DEVICE_PATH" ]; then
    echo "Error: baseline device template not found at ${BASE_DEVICE_PATH}." >&2
    exit 1
fi
if [ ! -f "$BASE_PULSES_PATH" ]; then
    echo "Error: baseline pulse template not found at ${BASE_PULSES_PATH}." >&2
    exit 1
fi

OUT_CSV="${REPO_DIR}/data/output/alg4/alg4_metrics.csv"
mkdir -p "$(dirname "$OUT_CSV")"
echo "algorithm,baseline_fidelity,merged_fidelity,baseline_latency_s,merged_latency_s,baseline_pulses,merged_pulses" >"$OUT_CSV"

for algo in qft shor adder w bv ghz; do
    if ! name="$(get_qasm_stub "$algo")"; then
        echo "Warning: unknown algorithm key '${algo}', skipping." >&2
        continue
    fi

    source_qasm="${REPO_DIR}/data/test_benchmark/experiments/${name}.qasm"
    if [ ! -f "$source_qasm" ]; then
        echo "Skipping ${algo}: source QASM ${source_qasm} missing." >&2
        continue
    fi

    alg_dir="${REPO_DIR}/data/output/alg4/${algo}"
    mkdir -p "$alg_dir"

    baseline_qasm="${alg_dir}/${name}_ideal_baseline.qasm"
    baseline_pulses="${alg_dir}/${name}_ideal_baseline_pulses.json"
    merged_template="${alg_dir}/${name}_ideal_merged_pulses.json"
    merged_device="${alg_dir}/${name}_ideal_merged_device.json"
    postmerge_qasm="${alg_dir}/${name}_ideal_postmerge.qasm"
    postmerge_pulses="${postmerge_qasm%.qasm}_pulses.json"
    baseline_sim="${alg_dir}/${name}_ideal_baseline_sim.json"
    postmerge_sim="${alg_dir}/${name}_ideal_postmerge_sim.json"

    echo "=== ${algo} (${name}) ==="

    # Step 1: compile baseline schedule from the logical circuit.
    rm -f "$baseline_qasm" "$baseline_pulses"
    ./build/QASMTrans \
        -i "$source_qasm" \
        -m rigetti \
        -c "$BASE_DEVICE_PATH" \
        -p "$BASE_PULSES_PATH" \
        -o "$baseline_qasm" \
        -v 0

    if [ ! -f "$baseline_pulses" ]; then
        echo "Error: baseline pulses not generated for ${algo}." >&2
        exit 1
    fi

    baseline_device=$(jq -r '.backend.backend_config // empty' "$baseline_pulses")
    if [ -z "$baseline_device" ] || [ "$baseline_device" = "null" ]; then
        baseline_device="$BASE_DEVICE_INPUT"
    fi
    baseline_device="$(resolve_path "$baseline_device")"

    # Step 2: optimise merged waveforms.
    rm -f "$merged_template" "$merged_device"
    echo "    optimising merged pulses (limit=${MERGE_LIMIT_INPUT})"
    python device_pulse_fab/merge_pulse.py \
        --pulse-dump "$baseline_pulses" \
        --pulse-template "$BASE_PULSES_PATH" \
        --device-config "$BASE_DEVICE_PATH" \
        --merge-limit "$MERGE_LIMIT_INPUT" \
        --output-pulses "$merged_template" \
        --output-device "$merged_device"

    if [ ! -f "$merged_template" ] || [ ! -f "$merged_device" ]; then
        echo "Error: merged artefacts not generated for ${algo}." >&2
        exit 1
    fi

    # Step 3: compile the post-merge schedule.
    rm -f "$postmerge_qasm" "$postmerge_pulses"
    ./build/QASMTrans \
        -i "$source_qasm" \
        -m rigetti \
        -c "$merged_device" \
        -p "$merged_template" \
        -o "$postmerge_qasm" \
        -v 0

    # Step 4: simulate baseline and merged schedules.
    rm -f "$baseline_sim"
    python device_pulse_fab/simulate_ankaa9q.py \
        --pulse "$baseline_pulses" \
        --device "$baseline_device" \
        --qasm "$baseline_qasm" \
        --output "$baseline_sim"

    rm -f "$postmerge_sim"
    python device_pulse_fab/simulate_ankaa9q.py \
        --pulse "$postmerge_pulses" \
        --device "$merged_device" \
        --qasm "$postmerge_qasm" \
        --output "$postmerge_sim"

    baseline_fid=$(jq -r '.fidelity' "$baseline_sim")
    merged_fid=$(jq -r '.fidelity' "$postmerge_sim")

    baseline_lat=$(jq -r '.backend.total_duration // (.schedule[-1].start_time + .schedule[-1].duration)' "$baseline_pulses")
    merged_lat=$(jq -r '.backend.total_duration // (.schedule[-1].start_time + .schedule[-1].duration)' "$postmerge_pulses")

    baseline_pulse_count=$(jq -r '.backend.total_pulses // (.schedule | length)' "$baseline_pulses")
    merged_pulse_count=$(jq -r '.backend.total_pulses // (.schedule | length)' "$postmerge_pulses")

    echo "${algo},${baseline_fid},${merged_fid},${baseline_lat},${merged_lat},${baseline_pulse_count},${merged_pulse_count}" >>"$OUT_CSV"
done

echo "Metrics written to ${OUT_CSV}"

if [[ $DO_PLOT -eq 1 ]]; then
    echo "Generating plots via scripts/plot_alg4_metrics.py"
    python "${REPO_DIR}/scripts/plot_alg4_metrics.py"
fi
