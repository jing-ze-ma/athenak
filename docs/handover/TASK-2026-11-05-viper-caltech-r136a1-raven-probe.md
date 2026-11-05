# TASK 2026-11-05 viper -> Caltech (USER go): send the R136a1 first-look setup, so viper can submit a 24 h probe on Raven NOW

USER: "submit the r136a1 probe now, ask caltech for it". viper has no R136a1 deck, IC or command line; none is on any fork branch.
Please push the files to this branch under `docs/handover/r136a1-1105/`. If a file is > ~50 MB, give its Caltech path instead and the user copies it.

## Please send
1. The deck: `.athinput` for the sp wedge, 512 x 128^2, 4 blocks, he_star_m1, fresh start.
2. The IC file(s) and tables the deck reads (EOS/opacity), with md5s.
3. The exact command line, with every command-line key you would launch with, and the output cadence.
4. The prod-1102 commit to build.
   - Default: viper uses the Raven binary already built for the gap runs, **23a55658 fmadoff (md5 620c8552)**.
   - Say if R136a1 needs a newer commit, e.g. 8cd7f0a9 (its defaults matter on fresh starts).
5. Anything known from the Caltech side: smoke results, s/cycle on H200, GPU memory per block, known failure modes.

## Plan on Raven (viper does this as soon as the files arrive)
- Layout: 1 node x 4 A100 = 4 ranks x 1 block. Estimate ~18 GB/block, scaled from RSG E: 19 GB for 548x64x64.
- Steps:
  1. a short smoke on 1 node, ~20-50 cycles;
  2. **one 24 h link** (Raven max wall is 24 h);
  3. stop rules as in the RSG/HG links;
  4. then report s/cycle, dt, sim time reached, and the NON-CONVERGED/floor counts.
- Run dir: `/raven/ptmp/jinma/raven_gap_1105/r136a1/` (ic/ and rst/ are ready).
- Fair share: Raven usage does not decay (resets monthly) and running jobs count at once. The probe is ~2 % of the RSG/HG gap usage.
  The aim is for it to finish before the RSG/HG chains queue on ~11-07, so the sooner the files arrive the better.

Reply with a NOTE on this branch.
