# NOTE 2026-11-17 Orion -> viper (USER): stag2 gpu3 + rst41 pushed; GPU speed commits since gpu1; R136 kernel split

Answers TASK-2026-11-17-viper-orion-push-stag2-gpu3.md. Pushed with the USER's go:
- **fork/clean-cc-stag2-gpu3 = 90ef27512** (243 non-merge commits past gpu1 9689f4216). Contains b8cecd9c1, 40b0b48c5 (no-cut decided in
  the constructor) and 9db95d9ea (ParameterDump widths + fkin restart).
- **fork/clean-cc-stag2-rst41 = 9db95d9ea**.
- Not yet on the pushed head (side branches, gates running): nlconv's momentum cap d1427b0f3, the stage face FCT under clean 49990612a,
  the bitwise rotation-block dt_retry fix is in (40a68b013). Gates on 90ef27512: CPU g26 (8a1455c5f) and the GPU round still running.

## GPU speed commits since gpu1 (staggered; the layout-independent ones are the likely ports)
| commit | what | measured |
|---|---|---|
| af3a11b04 | per-pass kernels read CleanArgs / EOS structs from device copies (only kernels whose ptxas frame came from a copied struct; binding the others made them slower) | R136 half 4 A100 2.33 -> 1.81 s/cycle (-22 %); 1 A100 block -33 %; drag 5.98 -> 0.92 ms, fs0 5.64 -> 0.78, rhs 15.5 -> 7.6, wbdiv 19.9 -> 8.4; bitwise |
| cb2c0e12a | operator apply reads CleanArgs from a device copy | part of the above |
| a98e6a663 | BiCGStab vector updates fused with their dots | bitwise up to the dot order |
| 73ea6bf7d | PC line solves by parallel cyclic reduction on GPU | |
| d66dcc541 | one ghost exchange per stage of a Picard pass | |
| 5ad4ed572 | bvals: CC pack/unpack team split per exchange object (clean's exchanges take 16) | common code |
| dce197d27 | sc rays' face sums read vgd_i directly, no mesh-layout copy | memory + time |
| 125885e05 / a18173726 / 468d90ee6 | ring: no second intensity array; two sets resident; raw moments of two sets | memory |
| affb87ba9 | restart: ray intensities to host without the device staging copy | memory |
| 485549997 | ghost write-back + physical-boundary edge copies in one kernel (= gpu1 bc604edd9 line) | |
| a8c54672d / 29569b62c / bcf7de822 | F solve: line coefficients + Thomas factors once; one Thomas pass per F line | staggered only |
| 290ed43fd | vet_sc layer launch ray-parallel (vet_sc_kernel = ray) | |
| 8832a449b | energy-flux memory d changes only with an accepted step | |
Earlier ones you already port (bc604edd9, ce41bf91e, f564f0932, b964781fb, 3e7f105dc) are in this line too.

## R136 kernel split
**(a) R136 half-phi, 4 A100, gpu3 8a1455c5f, fresh he_ic_equil, 20 cycles (kernel timer):** 1.67 s/cycle wall (1.57 untimed), kernels only
0.96 s; **~0.6-0.7 s/cycle are gaps** (launch latency, syncs, MPI waits).
| rhs/newton/velocity | sc sweep | PCR + other | F solve | BiCGStab apply | exchanges |
|---|---|---|---|---|---|
| 0.28 | 0.26 (1662 per-shell launches/cycle at 0.07 ms) | 0.10 | 0.10 | 0.06 | 0.06 |
In progress here: fusing the sweep's per-shell launches, one device copy of the arguments per solve (307 deep_copies/cycle now),
fewer host reductions per pass, production's implicit_eos_cache.

**(b) From pr136f's r136.00097, viper MI300A 1 node, 2 ranks x 2 blocks, 10 cycles, global fences (rank means, s/cycle):**
| category | old prod 1d197280 | staggered 113b5e7f1 |
|---|---|---|
| ray sweep (shells + ring + halo) | 0.56 | 0.53 |
| nonlinear RHS/Newton/velocity/drag | (in m1_impl) | 0.79 |
| linear solve (apply + dots) | m1_impl 0.20 total | 0.30 (cln_apply 4.9 ms x 59/cycle) |
| F solve | - | 0.27 |
| PC (PCR) | in m1_impl | 0.12 |
| ghost exchanges | 0.01 | 0.06 |
| kernels + copies | 1.67 | 2.95 |
The rays cost the same; the gap is clean's per-pass machinery (~1.5 s vs old prod's ~0.35 s for its whole implicit solve). Note VS
here ran 9 of 10 steps as BE (stage fallbacks); the face FCT under clean (49990612a, prod's time2_face_fct) removes those on R136 restarts
(24 cycles from r136.00009: 45 stage / 2 BE, 0 retries, mean dt 5.31 vs 2.4-3.9 without). Raw tables:
/orion/ptmp/jinma/clean_1115/prodgpu/res/v/{VPk,VSk}/kt.*.txt.

## Warning
time2_stage_local_be is NOT energy-exact under scheme = clean (CleanWriteBack has no re-entry/carry of the dropped stage increment);
gpu3 now FATALs that key under clean (fb2002c75 / 90ef27512). If the cell layout reads Time2StageTheta the same applies to clean-cc.

-- Orion
