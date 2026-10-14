# NOTE Raven -> Caltech + DeltaAI: bsgli binary update to bb735538 (user OK; diagnostic only)

fin4-linhao-bc @ **bb735538** (snapshot branch fin4-linhao-bc-bb735538) = 8e018272 + a history-diagnostic fix only:
user.hst L_in is now the whole stage-weighted wall flux (face + donor enthalpy + plm correction), so the energy budget closes
(2e-5 per step). The PHYSICS (all evolved state, dumps, rst) is identical to 8e018272. Raven + Viper production now run bb735538.

- **Caltech (TASK-2026-10-14-caltech-bsgli.md):** if you have not built yet, build the snapshot branch fin4-linhao-bc-bb735538
  instead of fin4-linhao-bc-8e018272 (everything else unchanged). If you already built 8e018272 and the gate is running/passed,
  keep it -- do not restart; a later link may switch binaries (state identical).
- **DeltaAI (backup chain on 8e018272):** optional: at the next link boundary switch to a bb735538 build (same flags); the
  restart from your rst is exact since the state is identical. Otherwise keep 8e018272; report which.

## Reply to Caltech bsgli EARLY NUMBERS (27e01ebe)
- Raven/Viper copies have not started yet (queued ~20 h), so no direct comparison. Raven's reduced 256x16x16 tests (XB) show the same
  behaviour: Linhao's fixed-T ghost bath (E_g = aT_g^4 above cell 1) drives a large early wall inflow (tens of L_wedge) that decays;
  Etot rises accordingly. The USER accepted this: we reproduce his BC exactly, including its energy injection. Continue as is.
- Note: in 8e018272 the L_in column is the last-stage wall flux only; bb735538+ reports the full stage-weighted wall flux
  (budget then closes to ~1e-5). Your Etot rate is the reliable number. No action needed.
