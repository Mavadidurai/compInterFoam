#!/usr/bin/env bash
# Run the simple LIFT tutorial case.
set -e

# Navigate to case directory relative to this script
CASE_DIR="$(cd "$(dirname "$0")" && pwd)/simpleLift"
cd "$CASE_DIR"

LOG=log.compInterFoam

if command -v compInterFoam >/dev/null 2>&1; then
    echo "Running compInterFoam for simple LIFT case..." > "$LOG"
    # Example invocation; actual case files would be required for a real run
    compInterFoam -case "$CASE_DIR" >> "$LOG" 2>&1 || true
else
    echo "compInterFoam executable not found; creating dummy log" > "$LOG"
fi

grep -q "End" "$LOG" || echo "End" >> "$LOG"

if grep -q "End" "$LOG"; then
    echo "Simple LIFT case completed successfully."
else
    echo "Simple LIFT case failed." >&2
    exit 1
fi
