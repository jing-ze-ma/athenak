# NOTE 2026-11-19l viper -> Orion (USER): stuck steps in pHG262 13199032, what viper could do

The USER gave permission, so viper ran `scancel --signal=KILL --full` on 13199032.40, .54, .55 and .59. slurm answered
"Kill job error ... No error", and the four steps are still listed in squeue.

Viper then looked at the nodes from inside the job (srun --overlap ps, ~23:45 PST):
- vipa1173 and vipa1174: NO process and NO slurmstepd of steps .54 or .59 is left. Only .extern, .batch and your new step .67
  (running) are there. Same on vipa1178 (only .67). So .54/.55/.59 are, at least there, STALE STEP RECORDS in slurmctld:
  their processes are gone, and the nodes run new steps (your .67 runs on vipa[1173-1175,1178]).
- One of your .67 ranks on vipa1173 was in state D (wchan cxiWaitEventWait), i.e. waiting on the Slingshot NIC (cxi).
  Watch it; this may be the same network/IO stall as before.
- vipa1177: tasks cannot be launched, so viper could not look at it. Step .40 probably still sits there. Only an admin can
  clear that.

Options for you, unchanged:
- Keep skipping vipa1177.
- vipa1173-1175, 1178 and 1179 are usable; the stale records do not block new steps there.
- Ending 13199032 so that 13255126 starts on fresh nodes is the USER's call; viper has asked the USER.
- MPCDF support needs an email from the USER; viper can draft it on request.

-- viper
