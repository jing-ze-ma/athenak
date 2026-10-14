# NOTE Raven -> Caltech + DeltaAI: bsgli binary update to bb735538 (user OK; diagnostic only)

fin4-linhao-bc @ **bb735538** (snapshot branch fin4-linhao-bc-bb735538) = 8e018272 + a history-diagnostic fix only:
user.hst L_in is now the whole stage-weighted wall flux (face + donor enthalpy + plm correction), so the energy budget closes
(2e-5 per step). The PHYSICS (all evolved state, dumps, rst) is identical to 8e018272. Raven + Viper production now run bb735538.

- **Caltech (TASK-2026-10-14-caltech-bsgli.md):** if you have not built yet, build the snapshot branch fin4-linhao-bc-bb735538
  instead of fin4-linhao-bc-8e018272 (everything else unchanged). If you already built 8e018272 and the gate is running/passed,
  keep it -- do not restart; a later link may switch binaries (state identical).
- **DeltaAI (backup chain on 8e018272):** optional: at the next link boundary switch to a bb735538 build (same flags); the
  restart from your rst is exact since the state is identical. Otherwise keep 8e018272; report which.
