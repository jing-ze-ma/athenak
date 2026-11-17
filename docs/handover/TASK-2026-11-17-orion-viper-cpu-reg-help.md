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

## UPDATE 1 (USER correction): use the GPU allocations, not CPU
USER: "I mean we have gpu allocations not cpu". So please run batch 1 on GPU, on the MI300A allocations you hold. A100 on Raven is optional,
since Orion also runs A100 itself. Keep it below your own gates, push 17 and production.
- Build each ref of the bundle (trimtest e0a52c637, clean-cc-stag2-gpu2-m36 5155aa332, clean-cc-stag2-gpu2-m12 680738756, and the bases
  e9227103a, ce49d1c1d) with your MI300A gate recipe (NOTE 8077e78a6), plus `-D Athena_M1_VGD_NMAX=512`.
- Run your MI300A gate deck set (gates1..4_mi300a style) with `clean_layout = staggered` where the decks allow it, branch vs base on the
  same GPU type, and include 3x repeats of 2 decks (determinism) and the RT decks (rot1s, csCnr: N/2 + restart + N/2 bitwise).
- Expected:
  - GPU vs GPU branch-vs-base: bitwise where the branch does not act (as on CPU);
  - m36 changes sp3 with poles / theta-periodic seams / cs sbimp + extrap1;
  - m12 changes only analytic-closure or SC-ray decks.
- Outputs: /viper/ptmp2/jinma/orion_stag2/req1/out/. NOTE: NOTE-2026-11-17-viper-orion-gpu-batch1.md.
- The CPU request above is withdrawn.

## UPDATE 2: new m12 head
The m12 history was rewritten (M1 is now restricted to analytic closures and clean_rays = sc). Use **clean-cc-stag2-gpu2-m12 = 7e6cc0bd4**
(ce49d1c1d + 5a1c7e973 + 7e6cc0bd4), not 680738756. The bundle in req1/ has been replaced and also carries the bases ce49d1c1d and e9227103a
as refs.
