---
name: restart-rot-potential-bug
description: pgen.cpp's RESTART constructor never read problem/rot_potential, so every restarted run with rot_potential + WB applied the radial centrifugal force TWICE (0.4 % of g); found 2026-09-09 because restart-based ablation arms all "survived" with +0.5 % mass, -90 % radial KE, 2x dt
metadata:
  type: project
---

Found 2026-09-09 (viper). Symptom: restart arms from the lhllc rot-4.0 and hllc rot-15.11
dumps never reproduced the vertex dt collapse, on the old and the new binary alike, and their
histories left the original within 0.1 rot: mass +0.5 %, 1-KE -70..-90 %, dt 12 -> 15-22 s.
A CPU restart showed the cycle-0 state BIT-IDENTICAL to the source dump, so the FORCING
differed. Audit: src/pgen/pgen.cpp's restart ctor reads the hot_jupiter_param fields but
omitted rot_potential (default false), while UserProblem builds the potential from pin (with
the centrifugal term). SourceFunc then keeps the explicit radial centrifugal term
(rotpot_src false) on top of the WB potential that already carries it: Omega^2 r/g = 4.3e-3.
Fix: one line (commit on 2026-09-09 after ae300379; see git log "restart constructor").

**Affected:** every CHAINED run on the 09-07 defaults after its first restart: cs_mhd_prod2
(chain from rot ~7), cs_hyd_rs/hllc and ausmpup (chains from rot ~8 to their collapse at 15.2),
any sp run with rot_potential=true restarted. lhllc and ppmx collapsed INSIDE their lead jobs
(no restart) -> the collapse itself does not need this bug. All first-round restart ablation
arms (cs_ablate_r/{ctl..cache1,h_*,old_*}) are INVALID.

**How to apply:** any parameter read in ProblemGenerator's from-scratch ctor must also be read
in the restart ctor (pgen.cpp keeps two copies of the list -- diff them whenever a field is
added). Gate a restart by comparing a continuous run and a run restarted mid-way for bitwise
history identity (wb_cache_every=1 to remove the cache phase). See [[wb-restart-cache-bug]]
(orion's zero-background kick, same class), [[restart-output-dt-trap]],
[[cs-vertex-dt-collapse-0907-defaults]].
