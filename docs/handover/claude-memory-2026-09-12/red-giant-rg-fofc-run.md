---
name: red-giant-rg-fofc-run
description: "RG_fofc 09-11 ~15:30: red giant FROM SCRATCH at nghost=3 (FOFC+PLM needs 3; all old rsts have 2 baked in) with fofc=true, vceil=0, efloor_from_ekin+dfloor_keep_velocity ON (user: keep the floor consistency fixes), keep_temperature off; chain 196662->196663->196664; pin13_fofc; 1.8x per cycle vs V9f"
metadata:
  type: project
---
Run dir /orion/ptmp/jinma/Athenak/red_giant/RG_fofc (rg.athinput, sub.sh, sub2.sh generic
chain link picking the newest rst; _probe/ has the failed-restart probes + the first launch's log).
Binary /orion/ptmp/jinma/Athenak/red_giant/build_impl/athena_pin13_fofc = current tree (FOFC
steps 3-5 + RT guard), built in athenak/build_fofc_rg. Config = I8_vceil's except: nghost 3,
fofc=true, vceil=0 (user: "don't use velocity ceiling"), efloor_from_ekin=true and
dfloor_keep_velocity=true (floor consistency fixes, kept by agreement), keep_temperature false,
rst every 1e5. Risks flagged by the launch agent: opac_tmin=-1 (clamp off) from t=0 is untested
(I3 was a restart test; V9f ran from scratch with 3200); rad_gate tuned at 5e5 now active from 0.
First launch (196659, vceil 5e7, fixes off) ran 1000 cycles fine: dt 30.67 flat, 1.80x cost per
cycle vs V9f (ghosts 1.37x predicted), mass exact, no counters fired; cancelled and resubmitted.
Expect link 1 ~3e5, link 2 ~6e5, link 3 to 9e5 (V9f mean dt ~12 s over 0-5.7e5).
**Why:** the FOFC plan's purpose was the top-cell negative density / vacuum-cell runaway; this run
tests whether FOFC replaces the ceiling. Deaths to beat: 5.7e5 shell launch, 5.87e5, 6.13e5.
**How to apply:** check out.txt (cycle/dt), rg.log (fofc column; first nonzero row = first flag;
eos_dfloor/efloor), rg.hydro.hst; analysis _analysis_0911/i8.py with RUN -> RG_fofc.

## 09-12 00:55: merged binary pin14_merge (c90bd869, build_merge_rg2, MPI+OpenMP) timed in
## RG_fofc/_timing (job 196697): 0.486 s/cycle vs pin13 0.535-0.545 = ~10% faster (viper's seam-pack
## + conduction-split commits); bit-identical on the red-giant gate -> RG_fofc/bin/athena SWAPPED to
## pin14 for links 2-3 (196663/196664); link 1 (196662) still runs pin13. Link 1 stdout is in _probe/out.txt.196659 (file moved after open); links 2-3 write out.txt. Mixed provenance is
## harmless (identical) but note it if a bitwise A/B against RG_fofc is ever needed.
