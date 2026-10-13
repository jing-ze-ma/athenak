# e1caltech-1013: e1chain_hook.sh = Raven jobs/e1chain_raven.sh (md5 bc4b5210, PS prodswitch-1013; validated BITWISE on A100 and
# MI300A) ported for Caltech. Only changes: E1C_W defaults to this script's directory (rst_pick.py, prune_rst.sh next to it), and an
# optional production gate (E1C_GATE=1): the link refuses to run unless $D/GO_e1 (copy of docs/handover/GO-e1-caltech.txt) AND
# $E1C_VERIFIED (written by smoke_gate.sh) exist. A refused link writes $D/NOGO and exits 1 (no STOP, no srun, no prune).
#
# CONTRACT (SOURCED by your link script right after it sources run.cfg; cwd = run dir):
#   in : D (run dir), IN (deck), KEYS (extra keys), TLIM, STOP_NC (default 1); optional E1C_NCYC (tests only), E1C_GATE, E1C_VERIFIED
#   out: RST (empty = fresh -i $IN; else -r $RST), KEYS += outputN/last_time=... time/tlim=$TLIM [time/nlim=...]
#   the link script then runs: srun ... $BIN ${RST:+-r $RST} (or -i $IN) -d $D -t <limit-10min> $KEYS > $D/run.log.$SLURM_JOB_ID; rc=$?
#   (rc MUST be a shell variable named rc; the EXIT trap below reads it; early exits never set rc and never prune)
# prodswitch-1013 (PS, USER GO 2026-10-13) e1chain_raven.sh: SOURCED by an E1 run.cfg inside a rc_1009/jobs/rc_raven.sh link
# (cwd = run dir D; rc_raven.sh then runs: srun $BIN ${RST:+-r $RST | -i $IN} -d $D -t <limit-10min> $KEYS; rc=$?).
# rc_raven.sh restarts only from RST= in run.cfg; this hook makes the E1 links chain like the BSG links (bsgrec_1011/jobs/link.sh):
#   DONE -> exit 0; STOP -> exit 1; GUARD: an older RUNNING job of the same name -> exit 0.
#   RST = the newest VALID rst in D/rst (jobs/rst_pick.py; fresh start from IN only if D/rst has no rst; rst files but none valid
#   -> STOP), every outputN/last_time = floor(t/dt)*dt (output cadence kept, as rst_info.py), time/tlim=$TLIM; rst t >= TLIM -> DONE.
#   E1C_NCYC (tests only): stop E1C_NCYC cycles after the start point.
#   EXIT trap (replaces the plain pruning trap): after a link whose srun ran: rc != 0, FATAL, NaN (last 60 lines), DIVERGED or
#   (STOP_NC=1) NON-CONVERGED > 0 -> STOP (no garbage restart chained); then prune_rst.sh "$D" 3 (keeps the newest 3; it
#   does nothing after a failed/STOPped link); newest valid rst t >= TLIM -> DONE.
# Needs in run.cfg: TLIM, STOP_NC (default 1), D (set by rc_raven.sh).
E1C_W=${E1C_W:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
E1C_PY=python3
: "${TLIM:?run.cfg must set TLIM}"; STOP_NC=${STOP_NC:-1}
D=${D:-$PWD}
[ -f "$D/DONE" ] && { echo "DONE: exit 0"; exit 0; }
[ -f "$D/STOP" ] && { echo "STOP ($(cat "$D/STOP")): exit 1"; exit 1; }
if [ "${E1C_GATE:-0}" = 1 ]; then
  if [ ! -f "$D/GO_e1" ] || [ ! -f "${E1C_VERIFIED:-/nonexistent}" ]; then
    echo "NOGO: GO_e1=$( [ -f "$D/GO_e1" ] && echo yes || echo NO ) VERIFIED=$( [ -f "${E1C_VERIFIED:-/nonexistent}" ] && echo yes || echo NO ): exit 1"
    echo "job $SLURM_JOB_ID refused: GO_e1 or smoke VERIFIED missing" >> "$D/NOGO"; exit 1
  fi
fi
if [ -n "$SLURM_JOB_ID" ]; then
  E1C_ME=$(squeue -h -j "$SLURM_JOB_ID" -o %S)
  for l in $(squeue -h -u "$USER" -t R -n "$SLURM_JOB_NAME" -o "%i,%S"); do
    id=${l%%,*}; st=${l#*,}; [ "$id" = "$SLURM_JOB_ID" ] && continue
    if [[ "$st" < "$E1C_ME" ]] || { [ "$st" = "$E1C_ME" ] && [ "$id" -lt "$SLURM_JOB_ID" ]; }; then
      echo "GUARD: older running job $id: exit 0"; exit 0
    fi
  done
fi
mkdir -p "$D/rst"
E1C_PICK=$($E1C_PY $E1C_W/rst_pick.py "$D/rst"); E1C_RC=$?
[ $E1C_RC -ne 0 ] && { echo "rst_pick rc=$E1C_RC (rst files but none valid): STOP"; echo "no valid rst (rst_pick rc $E1C_RC) in job $SLURM_JOB_ID" > "$D/STOP"; exit 1; }
if [ -n "$E1C_PICK" ]; then
  read -r RST E1C_T E1C_NC E1C_TL0 <<< "$(echo "$E1C_PICK" | head -1)"
  awk -v t="$E1C_T" -v l="$TLIM" 'BEGIN{exit !(t >= l*(1-1e-12))}' && { echo "rst t=$E1C_T >= tlim: DONE"; touch "$D/DONE"; exit 0; }
  KEYS="$KEYS $(echo "$E1C_PICK" | tail -1) time/tlim=$TLIM"
  echo "e1chain: restart from $RST (t = $E1C_T, cycle $E1C_NC)"
else
  RST=; E1C_NC=0; KEYS="$KEYS time/tlim=$TLIM"
  echo "e1chain: FRESH start from $IN"
fi
[ -n "${E1C_NCYC:-}" ] && KEYS="$KEYS time/nlim=$(( E1C_NC + E1C_NCYC ))"
e1chain_post() {
  [ -n "${rc:-}" ] || return 0          # early exits (CANCEL, md5, GUARD, DONE, STOP) never ran srun
  local L="$D/run.log.$SLURM_JOB_ID" why=""
  [ -f "$L" ] || L="$D/run.log"
  [ "$rc" != 0 ] && why="srun rc=$rc"
  grep -q FATAL "$L" && why="$why FATAL"
  tail -n 60 "$L" | grep -qiw nan && why="$why NaN"
  grep -q DIVERGED "$L" && why="$why DIVERGED"
  local ncv=$(grep -o "NON-CONVERGED=[0-9.e+-]*" "$L" | tail -1 | cut -d= -f2)
  [ "$STOP_NC" = 1 ] && [ -n "$ncv" ] && awk -v x="$ncv" 'BEGIN{exit !(x+0 > 0)}' && why="$why NON-CONVERGED=$ncv"
  if [ -n "$why" ]; then echo "e1chain: $why: STOP"; echo "$why in job $SLURM_JOB_ID" > "$D/STOP"; fi
  bash $E1C_W/prune_rst.sh "$D" 3 "$rc" "$SLURM_JOB_ID"
  if [ -z "$why" ] && [ -z "${E1C_NCYC:-}" ]; then
    local t2=$($E1C_PY $E1C_W/rst_pick.py "$D/rst" 2>/dev/null | head -1 | cut -d' ' -f2)
    [ -n "$t2" ] && awk -v t="$t2" -v l="$TLIM" 'BEGIN{exit !(t >= l*(1-1e-12))}' && { echo "e1chain: rst t=$t2 >= tlim: DONE"; touch "$D/DONE"; }
  fi
  return 0
}
trap e1chain_post EXIT
