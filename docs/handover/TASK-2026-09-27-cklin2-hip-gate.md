# TASK for viper: HIP gate of ck-lin2-fma (merged from Caltech 09-27)

`rt-integration` merge **8be0ad67** brings in branch `ck-lin2-fma` (5 commits, all in
`src/utils/two_stream_rt.hpp`). It was gated on Caltech (CPU + CUDA H100/H200) only. **HIP is not verified.**
The ck-next merge 6a9ffb0d is also still waiting for its HIP gate; both can be gated in one go.

## What changed

| commit | change |
|---|---|
| 6bf26403 | cin/cout formed in ck_lin_build and rt_chain_ck_jlin; ck_coef stores ck_ci/ck_co only when a reader still needs them (`ckcst_`) |
| b8d45c7e | lP layer weights (slots 2-4) formed with `CkLayW` in lin1p and jlin instead of stored |
| 3581844c | e0 formed with `CkE0` (expm1) in lin1p, lin_build and jlin; ck_c0 stored only under `ckcst_` |
| 1fa7e920 | ck_beam_tau as one team per (block, beam group, column), kappa rho and face radii in team scratch |
| 7a573de9 | `CkLayW`: under `__CUDA_ARCH__` only, `dtc = __fma_rn(0.5*kru, dzu, __dmul_rn(0.5*kro, dz))` |

About 7a573de9: the base `ck_lin_build` compiled `dtc = dt_l + dt_u` as `DFMA(0.5*kru, dzu, dt_l)`. Once lin1p
and jlin form the weights themselves, nvcc contracted that sum differently and broke GPU-bitwise. Spelling out
the fma restored it. **On HIP the fix is not active**, so amdclang may contract the sum its own way.

## Caltech results

- Timing, H200, nx1 256, production keys, interleaved arms, median of 8-cycle windows:
  - 2 GPUs: 23.93 -> 21.56 ms/cycle (-9.9 %).
  - 1 GPU: 40.55 -> 36.37 ms/cycle (-10.3 %).
- CPU A/B gate: bitwise for every commit.
- CUDA (H100), base vs 7a573de9: rst and hst bitwise. 1 vs 2 ranks: rst bitwise.
- Noise gate for 1fa7e920 (before the fma fix, i.e. not GPU-bitwise): 0.65x / 0.98x the noise-pair mean.
  This is the fallback argument if HIP is not bitwise.
- Records: `/resnick/groups/carnegie_poc/jingze/cklin2` on Caltech.

## To do on viper (apudev, gpu:2 per node)

1. HIP A/B, parent of the merge (732a777f) vs 8be0ad67, dhj production keys. Compare rst/bin, not hst bytes.
   - **Bitwise:** done; record it.
   - **Not bitwise:** either
     - find the HIP contraction the same way (single-site variants: ck_lin_build vs lin1p vs jlin forming
       `CkLayW`; then `__builtin_fma` forms under `__HIP_DEVICE_COMPILE__`), or
     - accept on a noise gate (>= 2-3 1e-14 kick members, deep eint rms 1-100 bar, pass <= ~1.2x).
     Ask the user which.
2. HIP 2 vs 1 rank: rst bitwise.
3. Optional: MI300A timing, same binary family, interleaved. The kernel-shape change in ck_beam_tau
   (team scratch ~18.5 KB) was tuned on H200 only.
