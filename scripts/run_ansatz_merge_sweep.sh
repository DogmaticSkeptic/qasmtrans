#!/usr/bin/env bash
# Sweep merge limits for new 10-layer QAOA/VQE ansatz circuits.

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"

CONDA_ACTIVATE="/opt/homebrew/anaconda3/bin/activate"
ENV_NAME="QASMtrans"
if [[ ! -f "$CONDA_ACTIVATE" ]]; then
    echo "Conda activate script not found at $CONDA_ACTIVATE" >&2
    exit 1
fi

# shellcheck disable=SC1091
source "$CONDA_ACTIVATE" "$ENV_NAME"

export MPLCONFIGDIR="${REPO_DIR}/.mplconfig"
mkdir -p "$MPLCONFIGDIR"

DEVICE_TEMPLATE="${REPO_DIR}/data/output/qft4_rigetti_ideal_device.json"
PULSE_TEMPLATE="${REPO_DIR}/data/output/qft4_rigetti_ideal_pulses.json"

if [[ ! -f "$DEVICE_TEMPLATE" || ! -f "$PULSE_TEMPLATE" ]]; then
    echo "Required template files not found: ${DEVICE_TEMPLATE} or ${PULSE_TEMPLATE}" >&2
    exit 1
fi

declare -a CONFIGS=(
    "qaoa10:qaoa4_l10_ansatz"
    "vqe10:vqe4_l10_ansatz"
)

for entry in "${CONFIGS[@]}"; do
    IFS=: read -r out_label base_name <<<"$entry"
    algo_dir="${REPO_DIR}/data/output/ansatz_sweep/${out_label}"
    qasm_src="${REPO_DIR}/data/test_benchmark/experiments/${base_name}.qasm"

    if [[ ! -f "$qasm_src" ]]; then
        echo "QASM source not found: $qasm_src" >&2
        exit 1
    fi

    mkdir -p "$algo_dir"

    baseline_prefix="${algo_dir}/${base_name}_baseline"
    baseline_qasm="${baseline_prefix}.qasm"
    baseline_pulses="${baseline_prefix}_pulses.json"
    baseline_sim="${baseline_prefix}_sim.json"
    metrics_csv="${algo_dir}/${out_label}_merge_metrics.csv"

    echo "=== Processing ${out_label} ==="
    echo "Running baseline QASMTrans..."
    rm -f "$baseline_qasm" "$baseline_pulses"
    ./build/QASMTrans \
        -i "$qasm_src" \
        -m rigetti \
        -c "$DEVICE_TEMPLATE" \
        -p "$PULSE_TEMPLATE" \
        -o "$baseline_qasm" \
        -v 0

    echo "Simulating baseline pulses..."
    rm -f "$baseline_sim"
    python device_pulse_fab/simulate_ankaa9q.py \
        --pulse "$baseline_pulses" \
        --device "$DEVICE_TEMPLATE" \
        --qasm "$baseline_qasm" \
        --output "$baseline_sim"
    base_fid=$(jq -r '.fidelity' "$baseline_sim")
    base_lat=$(jq -r '.backend.total_duration // (.schedule[-1].start_time + .schedule[-1].duration)' "$baseline_pulses")
    base_pulses=$(jq -r '.backend.total_pulses // (.schedule | length)' "$baseline_pulses")
    printf "  baseline -> fidelity %.6f, latency %.6e s, pulses %s\n" "$base_fid" "$base_lat" "$base_pulses"

    echo "merge_limit,fidelity,latency_s,total_pulses" >"$metrics_csv"
    printf "0,%s,%s,%s\n" "$base_fid" "$base_lat" "$base_pulses" >>"$metrics_csv"

    for limit in {1..10}; do
        echo "-- merge limit ${limit}"
        merged_prefix="${algo_dir}/${base_name}_merged_L${limit}"
        merged_pulses="${merged_prefix}_pulses.json"
        merged_device="${merged_prefix}_device.json"

        python device_pulse_fab/merge_pulse.py \
            --pulse-dump "$baseline_pulses" \
            --pulse-template "$PULSE_TEMPLATE" \
            --device-config "$DEVICE_TEMPLATE" \
            --output-pulses "$merged_pulses" \
            --output-device "$merged_device" \
            --merge-limit "$limit"

        postmerge_prefix="${algo_dir}/${base_name}_postmerge_L${limit}"
        postmerge_qasm="${postmerge_prefix}.qasm"
        postmerge_pulses="${postmerge_prefix}_pulses.json"
        postmerge_sim="${postmerge_prefix}_sim.json"

        ./build/QASMTrans \
            -i "$qasm_src" \
            -m rigetti \
            -c "$merged_device" \
            -p "$merged_pulses" \
            -o "$postmerge_qasm" \
            -v 0

        python device_pulse_fab/simulate_ankaa9q.py \
            --pulse "$postmerge_pulses" \
            --device "$merged_device" \
            --qasm "$postmerge_qasm" \
            --output "$postmerge_sim"

        fid=$(jq -r '.fidelity' "$postmerge_sim")
        lat=$(jq -r '.backend.total_duration // (.schedule[-1].start_time + .schedule[-1].duration)' "$postmerge_pulses")
        pulses=$(jq -r '.backend.total_pulses // (.schedule | length)' "$postmerge_pulses")
        printf "  L=%02d -> fidelity %.6f, latency %.6e s, pulses %s\n" "$limit" "$fid" "$lat" "$pulses"
        printf "%d,%s,%s,%s\n" "$limit" "$fid" "$lat" "$pulses" >>"$metrics_csv"
    done

    echo "Metrics written to $metrics_csv"
done
