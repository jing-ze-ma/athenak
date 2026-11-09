# TASK 2026-11-09 viper -> ORION: raycheap11 -- R alone is 1st order: stage order of the exchange (CPU)

Answer with `docs/handover/NOTE-2026-11-xx-orion-raycheap11.md`. raycheap10: with the hydro skipped (dbg_hydro_skip_t2),
R alone is still p ~ 1 on Ar4, and the R toggles were inert.

## viper code reading: how the gas energy is staged in hesdirk2 (rad_m1_time2.cpp Time2FormStage)
- The stage increment t2inc carries ALL channels M1_T2_NK: E, the gas momentum (3), the gas ENERGY (M1_T2_EN) and the
  face fluxes: the gas energy is integrated with the same SDIRK weights as E (K1, K2 slopes; c1 = a21 dt in stage 1,
  c31 dt, c32 dt in stage 2), not only by the write-back.
- Stage 2 starts from the Heun average 0.5 (U^n + Y2) for the radiation (m1_t2_avg) and, via the Heun hydro stage 2
  (gam0 = gam1 = 0.5), for the gas: consistent. The step is U^n + dt (0.5 K1 + b2 K2 + g K3) with b2 = 0.5 - g:
  2nd order for non-stiff problems (sum b c = 1/2).
- But the hesdirk2 stage 1 has c2 = 1 with a21 = 1 - g, a22 = g: sum_j a2j c_j = g != c2^2/2 = 1/2, i.e. STAGE ORDER 1.
  With a stiff exchange (k dt >> 1, Ar4's clump: chi = 50 against dt 0.04-0.3) a stage-order-1 SDIRK reduces to order 1
  (classical order reduction; the same as the 09-26 note "coupling split 1st order when the exchange is stiff").
- The code already has `rad_m1/time2_tableau = trbdf2` (TR-BDF2: c = (0, 2g, 1), A = [[0,0,0],[g,g,0],[w,w,g]],
  STAGE ORDER 2; its explicit partner is the RK2 with c2 = 2 - sqrt 2, set by the driver). That is the candidate fix.

## Binary (REBUILD)
MOD = `fork/raycheap2-1109` head `addb642d0` (test forcing only: the stage-1 abscissa of trbdf2, c2 = 2g, for the
modulated inner flux; dbg_hydro_skip_t2 as in raycheap10). viper compile-checked it (incremental MI300A build).

## Runs (your Ar4 deck, t_end 10, 4 levels cfl 0.15 / 0.075 / 0.0375 / 0.01875, 24 ranks; all with
## rad_m1/dbg_hydro_skip_t2=true unless noted)
| arm | overrides | question |
|---|---|---|
| F1 | time2_vet_col=lag vet_col_every=1000000 | closure fully frozen: exchange + transport staging alone |
| E1 | closure=eddington implicit_hr_model=true | the no-VET path |
| N1 | implicit_opac_update=true implicit_gas_newton=true implicit_opac_newton=true time2_vet_col=predict | "Newton opacity" (B2 done properly) |
| T1 | coupling=false time2_vet_col=predict | no exchange: transport alone with the time-varying boundary (expect p = 2) |
| R1 | time2_tableau=trbdf2 time2_vet_col=lag vet_col_every=1000000 | stage order 2 + frozen closure (expect p = 2 if the stage order is the term) |
| R2 | time2_tableau=trbdf2 | stage order 2 with the deck's update (extrap2) |
| R3 | time2_tableau=trbdf2, WITHOUT dbg_hydro_skip_t2 | the full coupled problem with TR-BDF2 |
(The deck's keys use the rad_m1/ prefix on the command line. If a key combination is refused, report and skip.)

## Analysis
Per arm: conv.py d, p (finest pair), walls, hesdirk2 / BE / fallback counts. Key reading: T1 ~2 confirms transport;
F1 ~1 with R1 ~2 confirms the stage-order reduction of the exchange as THE term; R3 tells whether TR-BDF2 fixes the
full problem.

## NOTE back
Per-arm conv.py lines with counts and walls. No code changes.
