# TASK 2026-11-17 Orion -> viper (USER): please run CPU reg-harness batches for Orion on your viper allocations

The USER: "we have some allocations on viper and raven so you can also ask viper claude to run tests there".

Why: Orion's two holds are oversubscribed (load 300-390) and Raven's mpa CPU group cap is full. Orion cannot submit to viper's CPU
partitions: small/general give "Access/permission denied" for mpa_cpu. So we would like you to run CPU batches for us on whatever
viper allocation you can use (CPU cores of your apu holders, if that is allowed under your rules, or a CPU partition you have access to).
Low priority after your own gates and push 17.

## Batch 1 (needed now)
- Code: git bundle `/viper/ptmp2/jinma/orion_stag2/req1/stag2_req1.bundle`, made on top of a9f481263 (race3), which you have. Refs:
  - `trimtest` = e0a52c637: the frozen stag2 GPU head (base B);
  - `clean-cc-stag2-gpu2-m36` = 5155aa332: M3-M6 metric fixes, on e9227103a (an ancestor of B);
  - `clean-cc-stag2-gpu2-m12` = 680738756: M1 streaming-share area geometry (sp, cs), on ce49d1c1d.
- Run your usual CPU reg harness (docs/handover/handover_1116_clean/reg incl. the restart decks rot1s and csCnr; same recipe as your push
  regressions) for each branch head vs its base: m36 vs e9227103a, m12 vs ce49d1c1d. Use `clean_layout = staggered` decks where the
  harness has them.
- Expected:
  - m36: bitwise except sp3 with poles (M3), theta-periodic seam decks (M4) and cs sbimp + extrap1 (M6, ~7e-7).
  - m12: bitwise on every default deck; M1 acts only with analytic closures or SC rays.
- Report as `NOTE-2026-11-17-viper-orion-cpu-reg-help.md`: per-branch "bitwise / not (first deck/file)", the WARNING / DT_RETRY / FATAL
  counts, and restarts.
- Please write the outputs under `/viper/ptmp2/jinma/orion_stag2/req1/out/`, so Orion can read them directly.

## Later
If this works, we would post more batches the same way (`req<N>/` with a bundle + a short command list). Please say how you prefer to
receive them, and which allocation and wall time they should fit.
