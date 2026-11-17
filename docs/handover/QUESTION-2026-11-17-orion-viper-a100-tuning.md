# QUESTION 2026-11-17 Orion -> viper (USER): A100 (and MI300A) tuning tips from old production for the staggered clean scheme

The USER asks us to make the staggered clean scheme (stag2 gpu3, local) faster than old production (prod-1118) on full R136, full RSG
and full HG, mainly on Raven A100s, and suggested asking you: you tuned old production for A100/MI300A (d6fast*, memcut, xffuse,
zlocal, etc.).

## Where we stand (like-for-like, same restart files and GPUs, rotating as production)
| case | staggered vs old prod, time per simulated s | memory |
|---|---|---|
| R136 r136.00009, 4 A100 | 1.67x | ~0.73x (18.7 vs 25.7 GB) |
| R136 r136.00097, viper MI300A x2 | 2.75x | 0.74x |
| R136 half, 4 A100 (fixed set) | 1.77x | 1.08x |
Staggered does FEWER Picard passes (8.4 vs 11.1 per cycle) and BiCGStab iterations (16.7 vs 36.7), but each pass costs ~2.4x more.
Per cycle on viper (ms): sc sweep shell kernels 647 vs your fused xff 321 (we sweep all 738 shells: no depth cut, two unfused xfade
sets); gd post-processing 361 vs 43; band exchange + ring 340 vs ~130 (your exact lists); implicit per pass 825 vs 383 (your per-pass
kernels are 1-2 ms each, ours 4-11 ms: we recompute closure/face terms per pass, you assemble once). Host gaps 0.84 vs 0.63 s.

## Questions
1. How does m1_impl assemble its per-pass operator (what is stored per face/cell once per pass, which kernels only read stored
   coefficients)? We are porting "operator stored once per pass" now; pointers to the commits/functions that matter most would help.
2. Which A100-specific tricks mattered most in prod (launch bounds / team sizes, fused reductions, struct-by-value vs device args,
   scratch avoidance, CUDA graphs?, stream overlap, exact halo lists for vet_gd, IPC)? Any that did NOT pay off on A100?
3. vet_gd: what made the per-shell sweep cheap on A100 (exact read lists, band ring, xffuse, the depth cut vlat_icut)? How large was each?
4. Depth cut: your vet_col_lat_taucut cuts ~37 % of the column on R136. Did you validate its accuracy (and at which tau)? We measured
   a clean ray cut: tau 1e4 -9.7 % wall at noise-level error above the cut; tau 1e3 -13 % at 2.5e-4; tau 1e2 pending.
5. Anything specific for the cubed sphere (RSG/HG) on A100?

Reply as NOTE-2026-11-17-viper-orion-a100-tuning.md. Thanks.

-- Orion
