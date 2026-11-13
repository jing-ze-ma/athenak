# NOTE 2026-11-13 15:55 PST Caltech -> viper: T9L6.2 UPDATE3 running on 2 H200 (short2g); restart NOT bitwise; dhj 30 G status

## Binary

- prod-1118 1d1972805, he_star_m1, CUDA Release, MPI, `-ffp-contract=off` + fmad off, HOPPER90. Incremental build of
  the p1103 tree: 235 files recompiled.
- `builds_rstfixx/bin/athena_gpu_he_star_m1_prod1118_1d197280_fmadoff`, md5 26a473d92edf5e38c8c2626d44db1cce.

## Run setup (t9l62_ourwall.athinput only, UPDATE3 overrides)

- Dir: `/resnick/groups/carnegie_poc/jingze/t9l62_1113`.
  - Inputs in `in/`: md5_unpacked.txt checks OK on all 4 files.
  - T9DIR was replaced in the deck.
- Layout: 512 x 128 x 256, theta 3pi/8..5pi/8, phi 0..pi/2, 2 MeshBlocks.
- Placement: 2 ranks x 1 H200 in short2g (job 6854764, hpc-sm-02-03).

**Wall L_in: no scaling applied, none needed.**
- `<rad_m1>/implicit_flux_x1min` is a FLUX (he_star_m1.cpp:32: `implicit_flux_x1min = L/(4 pi r_in^2)`), so it is
  independent of the wedge.
- Checked in the run: the history's L_in (wedge luminosity, F_in A(r_in)) = 2.371e38 erg/s
  = 2.480e39 (L_IC = 6.9466e13 x 4 pi r_in^2) x 0.0957 (wedge solid angle 1.202 sr / 4 pi).

## Smoke (concurrent with the production run on the same 2 GPUs)

- 10 cycles: rc 0. 20 cycles fresh: rc 0. 10-cycle restart from the cycle-10 rst: rc 0, no FATAL.
- **Restart is NOT bitwise.** Checked with bin_double=true in a smoke copy of the deck (`t9l62_smoke.athinput`;
  `bin_double` cannot be added on the command line).
  - First divergence: cycle 11, the first step after the restart. The DT_RETRY residual is 3.642358e-04 (fresh) vs
    3.68999e-04 (restart).
  - After that the dt histories split. Fresh: one retry at cycle 11. Restart: retries at cycles 11, 16 and 19.
  - At cycle 20, t = 2011.68 s (fresh) vs 1865.68 s (restart). Every field differs (dens rel 3.6e-4, m1_f1 O(1)).
  - The restart log prints `rstabe-1109: time2_auto_be skip run restored (len 1, left 0)`. The cycle-11 residual
    already differs, so some M1 solver state that affects the first post-restart solve is not in the rst, or is not
    restored. Candidates: lag / Picard / VET cache, or the time2_posfix / auto_be state around the cycle-8 fallback.
  - Not tested: fresh-vs-fresh determinism on GPU for this deck, so GPU non-determinism is not excluded. A CPU
    np-1/np-2 restart gate would settle it. No such test was run here.
  - Production is not affected while it runs straight through. But every link restart will NOT be bitwise until
    this is understood.

## Production (UPDATE3: run to slot end - 10 min, rst at the end)

- Fresh start 15:39 PST, `-t 1:00:56`, so it stops ~16:40 PST and writes its final rst in `prod/rst/`.
- Early warnings, all recovered:
  - cycle 0: 9.4 M of 16.8 M cells hit eos_dfloor, presumably the low-density outer IC (not checked);
  - cycles 8, 13, 16: hesdirk2 stage 1 NOT POSITIVE -> time2_posfix, stage 2 NOT ADMISSIBLE -> BE redo;
  - DT_RETRY at cycles 11 and 64 (NON-CONVERGED, one cut each); the cap was released at cycle 79.
- dt = 129.9-130.0 s (CFL-limited) after the start-up.
- Speed: cycle 120 reached at t = 4.0 h simulated, 7.6 min wall after the start. That is ~3.8 s/cycle including init,
  with the GPUs shared with the smokes for ~6 of those minutes. A clean number will be in the final log.
  - At 130 s/cycle and ~3.5 s/cycle, 800 h = ~22,000 cycles = ~21 h on 2 H200 for this wedge.

L(r_out) (wedge, erg/s; L_top = outermost shell) and the mass flux at the top:

| t (h) | L_bot | L_top | L_in | Mdot_top (g/s, wedge) |
|---|---|---|---|---|
| 1.0 | 2.227e38 | 2.349e38 | 2.371e38 | 4.55e14 |
| 3.0 | 2.329e38 | 2.353e38 | 2.371e38 | 4.37e14 |
| 5.0 | 2.369e38 | 2.360e38 | 2.371e38 | 4.18e14 |

L_top / L_in = 0.995 at 5 h. The final rst and a longer L(t) will follow after the slot ends (~16:40 PST).

## Next slot

- short2g ends 16:55 PST.
- The continuation from the final rst goes into the next 2-GPU slot: test2g 6854659 (3 d) is pending, no estimate.
- No 4-GPU short holder was submitted (per UPDATE3).

## dhj (TOI-1518b 20 G only; WASP-189b 30 G cancelled for Caltech)

- The dhj prod-1118 build job 6862110 (expansion, CPU) is pending in the queue. Low priority.
- The bundle caltech_xfer_1113 has NOT arrived at Caltech. Nothing dhj has started.
