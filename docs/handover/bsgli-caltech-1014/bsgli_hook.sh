# bsgli-caltech-1014 bsgli_hook.sh: the chain hook of the Caltech Linhao-IC BSG production (bsgLIc).
# = e1caltech-1013/e1chain_hook.sh (md5 3fb8e2dd; Raven prodswitch-1013 e1chain_raven.sh, validated BITWISE on A100 and MI300A, and
# run by your E1 chain) with the variables renamed E1C_ -> BLC_ and ONE change: the production gate needs only the smoke VERIFIED
# file (the USER approved this production: NO GO file). BLC_GATE=1 and no $BLC_VERIFIED -> $D/NOGO, exit 1 (no srun, no STOP, no prune).
#
# CONTRACT (SOURCED by link_bsgli.sh right after it sources run.cfg; cwd = run dir):
#   in : D (run dir), IN (deck), KEYS (extra keys), TLIM, STOP_NC (default 1); optional BLC_NCYC (tests only), BLC_GATE, BLC_VERIFIED
#   out: RST (empty = fresh -i $IN; else -r $RST), KEYS += outputN/last_time=... time/tlim=$TLIM [time/nlim=...]
#   (rc MUST be a shell variable named rc; the EXIT trap reads it; early exits never set rc and never prune)
# Rules: DONE -> exit 0; STOP -> exit 1; GUARD: an older RUNNING job of the same name -> exit 0.
#   RST = the newest VALID rst in D/rst (rst_pick.py; fresh start from IN only if D/rst has no rst; rst files but none valid -> STOP),
#   every outputN/last_time = floor(t/dt)*dt (output cadence kept), time/tlim=$TLIM; rst t >= TLIM -> DONE.
#   EXIT trap after a link whose srun ran: rc != 0, FATAL, NaN (last 60 lines), DIVERGED or (STOP_NC=1) NON-CONVERGED > 0 -> STOP;
#   then prune_rst.sh "$D" 3 (keeps the newest 3; does nothing after a failed/STOPped link); newest valid rst t >= TLIM -> DONE.
BLC_W=${BLC_W:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
BLC_PY=python3
: "${TLIM:?run.cfg must set TLIM}"; STOP_NC=${STOP_NC:-1}
D=${D:-$PWD}
[ -f "$D/DONE" ] && { echo "DONE: exit 0"; exit 0; }
[ -f "$D/STOP" ] && { echo "STOP ($(cat "$D/STOP")): exit 1"; exit 1; }
if [ "${BLC_GATE:-0}" = 1 ] && [ ! -f "${BLC_VERIFIED:-/nonexistent}" ]; then
  echo "NOGO: smoke VERIFIED missing (${BLC_VERIFIED:-unset}): exit 1"
  echo "job $SLURM_JOB_ID refused: smoke VERIFIED missing" >> "$D/NOGO"; exit 1
fi
if [ -n "$SLURM_JOB_ID" ]; then
  BLC_ME=$(squeue -h -j "$SLURM_JOB_ID" -o %S)
  for l in $(squeue -h -u "$USER" -t R -n "$SLURM_JOB_NAME" -o "%i,%S"); do
    id=${l%%,*}; st=${l#*,}; [ "$id" = "$SLURM_JOB_ID" ] && continue
    if [[ "$st" < "$BLC_ME" ]] || { [ "$st" = "$BLC_ME" ] && [ "$id" -lt "$SLURM_JOB_ID" ]; }; then
      echo "GUARD: older running job $id: exit 0"; exit 0
    fi
  done
fi
mkdir -p "$D/rst"
BLC_PICK=$($BLC_PY $BLC_W/rst_pick.py "$D/rst"); BLC_RC=$?
[ $BLC_RC -ne 0 ] && { echo "rst_pick rc=$BLC_RC (rst files but none valid): STOP"; echo "no valid rst (rst_pick rc $BLC_RC) in job $SLURM_JOB_ID" > "$D/STOP"; exit 1; }
if [ -n "$BLC_PICK" ]; then
  read -r RST BLC_T BLC_NC BLC_TL0 <<< "$(echo "$BLC_PICK" | head -1)"
  awk -v t="$BLC_T" -v l="$TLIM" 'BEGIN{exit !(t >= l*(1-1e-12))}' && { echo "rst t=$BLC_T >= tlim: DONE"; touch "$D/DONE"; exit 0; }
  KEYS="$KEYS $(echo "$BLC_PICK" | tail -1) time/tlim=$TLIM"
  echo "bsglichain: restart from $RST (t = $BLC_T, cycle $BLC_NC)"
else
  RST=; BLC_NC=0; KEYS="$KEYS time/tlim=$TLIM"
  echo "bsglichain: FRESH start from $IN"
fi
[ -n "${BLC_NCYC:-}" ] && KEYS="$KEYS time/nlim=$(( BLC_NC + BLC_NCYC ))"
bsglichain_post() {
  [ -n "${rc:-}" ] || return 0          # early exits (CANCEL, md5, GUARD, DONE, STOP) never ran srun
  local L="$D/run.log.$SLURM_JOB_ID" why=""
  [ -f "$L" ] || L="$D/run.log"
  [ "$rc" != 0 ] && why="srun rc=$rc"
  grep -q FATAL "$L" && why="$why FATAL"
  tail -n 60 "$L" | grep -qiw nan && why="$why NaN"
  grep -q DIVERGED "$L" && why="$why DIVERGED"
  local ncv=$(grep -o "NON-CONVERGED=[0-9.e+-]*" "$L" | tail -1 | cut -d= -f2)
  [ "$STOP_NC" = 1 ] && [ -n "$ncv" ] && awk -v x="$ncv" 'BEGIN{exit !(x+0 > 0)}' && why="$why NON-CONVERGED=$ncv"
  if [ -n "$why" ]; then echo "bsglichain: $why: STOP"; echo "$why in job $SLURM_JOB_ID" > "$D/STOP"; fi
  bash $BLC_W/prune_rst.sh "$D" 3 "$rc" "$SLURM_JOB_ID"
  if [ -z "$why" ] && [ -z "${BLC_NCYC:-}" ]; then
    local t2=$($BLC_PY $BLC_W/rst_pick.py "$D/rst" 2>/dev/null | head -1 | cut -d' ' -f2)
    [ -n "$t2" ] && awk -v t="$t2" -v l="$TLIM" 'BEGIN{exit !(t >= l*(1-1e-12))}' && { echo "bsglichain: rst t=$t2 >= tlim: DONE"; touch "$D/DONE"; }
  fi
  return 0
}
trap bsglichain_post EXIT
