#!/bin/bash
# Run an axisymmetric fs-LIFT case in three stages, as used for the reference result:
#   stage 1:  0 - 2 ps      time step at most 10 fs  (laser pulse and electron-lattice coupling)
#   stage 2:  2 ps - 1 ns   time step at most 1 ps   (film launch), fields saved every 20 ps
#   stage 3:  1 ns - END    time step at most 2.5 ps (flight), fields saved every 0.2 ns
# The time step is otherwise set by the Courant limits in system/controlDict.
#
# Usage:  tools/run_case.sh cases/<case> [number of processes, default 4] [end time in s, default 2e-8]
#         RESUME=1 tools/run_case.sh cases/<case> ...   continue from the latest saved time
# Droplet detection on the beam axis is only reported (run_status.txt, droplet_detections.txt);
# it never stops the run. Confirm any detection in ParaView.
# Optional environment variables (defaults are the reference values above):
#   MAXDT2, MAXDT3    largest time step of stages 2 and 3 [s] (a finer mesh needs smaller ones,
#                     see tools/make_mesh_level.py)
#   WRITE1, WRITE2, WRITE3   write interval of stages 1, 2 and 3 [s]
set -e
TOOLS="$(cd "$(dirname "$0")" && pwd)"
CASE="$(cd "${1:?usage: tools/run_case.sh cases/<case> [NP] [END]}" && pwd)"
NP=${2:-4}
END=${3:-2e-8}
MAXDT2=${MAXDT2:-1e-12}; MAXDT3=${MAXDT3:-2.5e-12}
WRITE2=${WRITE2:-2e-11}; WRITE3=${WRITE3:-2e-10}
cd "$CASE"
: "${FOAM_BASHRC:=/usr/lib/openfoam/openfoam2406/etc/bashrc}"
[ -n "$WM_PROJECT" ] || source "$FOAM_BASHRC"
if [ "$NP" = 1 ]; then RUN="compInterFoam"; else RUN="mpirun -np $NP compInterFoam -parallel"; fi
touch droplet_report_only

set_ctrl() {  # set_ctrl <endTime> <maxDeltaT> <writeInterval>
  sed -i -E "s/^endTime .*/endTime         $1;/; s/^maxDeltaT .*/maxDeltaT       $2;/; s/^writeInterval .*/writeInterval   $3;/" system/controlDict
}
latest() {  # latest saved time (processor0 for parallel runs)
  local d=.; [ "$NP" != 1 ] && d=processor0
  ls -d "$d"/[0-9]* 2>/dev/null | xargs -n1 basename | sort -g | tail -1
}
run_stage() {  # run_stage <name> ; runs in the background and reports droplets every minute
  echo "$1 start $(date -u +%T)" >> run_status.txt
  if [ -f "log.$1" ]; then  # keep the log of an interrupted stage
    local n=1; while [ -f "log.$1.part$n" ]; do n=$((n+1)); done; mv "log.$1" "log.$1.part$n"
  fi
  $RUN > "log.$1" 2>&1 &
  local pid=$!
  while kill -0 $pid 2>/dev/null; do
    sleep 60
    echo "$(date -u +%T) $(python3 "$TOOLS/detect_droplet.py" . 2>/dev/null)" >> run_status.txt
  done
  local rc=0; wait $pid || rc=$?
  echo "$1 exit $rc $(date -u +%T)" >> run_status.txt
  return $rc
}

if [ -z "${RESUME:-}" ]; then
  # mesh, initial fields and decomposition
  blockMesh > log.blockMesh 2>&1
  rm -rf 0 && cp -r 0.orig 0
  topoSet > log.topoSet 2>&1
  setFields > log.setFields 2>&1
  if [ "$NP" != 1 ]; then
    sed -i -E "s/^(numberOfSubdomains[[:space:]]+)[0-9]+;/\1$NP;/" system/decomposeParDict
    decomposePar -force > log.decomposePar 2>&1
  fi
  set_ctrl 2e-12 1e-14 "${WRITE1:-5e-13}";  run_stage stage1
  set_ctrl 1e-9  "$MAXDT2" "$WRITE2";  run_stage stage2
else
  t=$(latest); echo "resume from $t $(date -u +%T)" >> run_status.txt
  if python3 -c "import sys; sys.exit(0 if float('$t') < 1e-9 else 1)"; then
    set_ctrl 1e-9 "$MAXDT2" "$WRITE2"; run_stage stage2
  fi
fi
if python3 -c "import sys; sys.exit(0 if float('$END') > 1.0001e-9 else 1)"; then
  set_ctrl "$END" "$MAXDT3" "$WRITE3"; run_stage stage3
fi
