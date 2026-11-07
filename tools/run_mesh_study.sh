#!/bin/bash
# Mesh study of cases/fsLIFT_332mJ_Ti_meshStudy on three meshes (tools/make_mesh_level.py).
#
#   tools/run_mesh_study.sh <work dir> [number of processes, default 4]
#
# Part 1: levels 1, 2 and 3 from 0 to 1 ns (laser heating and launch of the film), one after the
#         other, each with the largest time steps scaled to its cell size.
# Part 2: levels 2 and 1 continued from 1 ns to 20 ns (flight and droplet formation).
# Level 3 is not continued past 1 ns: on 4 processes it would need about 45 hours more.
# Each level runs in <work dir>/L<level>; progress is in <work dir>/study_status.txt.
set -e
TOOLS="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$TOOLS/.." && pwd)"
WORK="$(mkdir -p "${1:?usage: tools/run_mesh_study.sh <work dir> [NP]}" && cd "$1" && pwd)"
NP=${2:-4}
STATUS="$WORK/study_status.txt"
say() { echo "$(date -u +%FT%T) $*" | tee -a "$STATUS"; }
disk_ok() {  # stop before the disk fills: a full disk corrupts the written fields
  local avail; avail=$(df --output=avail -m "$WORK" | tail -1)
  if [ "$avail" -lt "${MIN_FREE_MB:-500}" ]; then say "stop: only ${avail} MB free"; exit 3; fi
}
setup() {  # setup <level>
  local d="$WORK/L$1"
  if [ ! -d "$d" ]; then
    cp -r "$REPO/cases/fsLIFT_332mJ_Ti_meshStudy" "$d"
    python3 "$TOOLS/make_mesh_level.py" "$d" "$1" | tail -1 > "$d/maxdt.env"
  fi
}
for L in 1 2 3; do
  setup $L; disk_ok
  if grep -q "^End" "$WORK/L$L/log.stage2" 2>/dev/null; then say "L$L already at 1 ns"; continue; fi
  source "$WORK/L$L/maxdt.env"
  say "L$L part 1 start (MAXDT2=$MAXDT2)"
  MAXDT2=$MAXDT2 WRITE1=1e-12 WRITE2=2e-10 bash "$TOOLS/run_case.sh" "$WORK/L$L" "$NP" 1e-9
  say "L$L part 1 done"
done
for L in 2 1; do
  disk_ok
  source "$WORK/L$L/maxdt.env"
  say "L$L part 2 start (MAXDT3=$MAXDT3)"
  RESUME=1 MAXDT3=$MAXDT3 WRITE3=5e-10 bash "$TOOLS/run_case.sh" "$WORK/L$L" "$NP" 2e-8
  say "L$L part 2 done"
done
say "mesh study finished"
