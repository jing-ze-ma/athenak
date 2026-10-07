# Plaskett stream-impact page (DeltaAI)

Data stay on DeltaAI. DeltaAI makes the figures with `mkfigs_plaskett.py` and publishes the
page itself; only a NOTE line and the two numbers JSONs come back to `accretor-1006`.

## 1. Figures

Python 3 with numpy, scipy, matplotlib (any recent conda/module python; nothing else).
`<ak>` = your athenak checkout of `accretor-1006` (the script uses `<ak>/vis/python/bin_convert.py`).

```bash
cd <ak> || exit 1
git pull fork accretor-1006
python3 docs/handover/plaskett-1007/page/mkfigs_plaskett.py \
  --run <run>/q12_s1 <run>/q12f_s1 --vis <ak>/vis/python \
  --out <ak>/docs/handover/plaskett-1007/results
```

Output: `fig_q12_s1_maps.png`, `fig_q12_s1_budget.png`, the same for `q12f_s1`,
`fig_plaskett_cmp.png`, `numbers_q12_s1.json`, `numbers_q12f_s1.json`. Parameters come from
the input embedded in the first restart (`<run>/<arm>/rst/*.rst`, includes the command-line
keys), history columns by label from `ryper.user.hst`, FATAL count and min dt from every
`*.log`/`*.out`/`*.err` in the run dir (put the slurm log there). Look at the PNGs before
publishing. Tested on viper q12_s45 / q12_a1 / q12f_a1 (and a synthetic restart-appended hst).

## 2. Page (Artifact tool, new artifact)

Short HTML page, title **"Plaskett Stream Impact"**, the six PNGs as supporting files
(`files`), sections:
- **Setup** (4-6 lines, from `docs/handover/TASK-2026-10-07-plaskett-remote.md` and the
  inputs): Plaskett progenitor at the start of accretion, M_a 16 / M_d 18 Msun, a 32.56 Rsun,
  P 3.69 d (P_orb 0.4583 code), R_acc 9.0 Rsun; ry_per_accretor envelope mode (n = 3 envelope,
  c_ph 20 km/s, rho_ph 0.3 rho_stream, hot ambient 300 km/s); grid r 6.3-13.50 Rsun
  (0.85 d_L1), 500 x 4 x 2048, polynomial+plateau radial stretch, dr ~2.6e-3 Rsun at the
  surface; half an orbit. **env12 vs env12f**: same well-balancing (only below
  wb_rmax 8.52 = R_acc - 1.5 penetration depths), env12 = PLM + hllc, nghost 2;
  env12f = the same + FOFC (nghost 3).
- **Maps** per arm (`fig_<arm>_maps.png`): theta-mean density in units of the stream peak,
  full domain (top) and the impact zoom (bottom, cyan box); white dashed = ballistic L1
  orbit, grey disk = inner boundary, blue = R_acc.
- **Budget** per arm (`fig_<arm>_budget.png`).
- **Comparison** (`fig_plaskett_cmp.png`): accreted mass and cumulative specific AM through
  R_acc in j_Kep(R_acc), with the ballistic and star-rotation values.
- **Numbers**: a table from the two JSONs (t_final_orbits, M_stream_in_cum,
  M_through_Racc_cum, j_through_Racc_last0p1orb_over_jK, j_ballistic_over_jK, fatal_count,
  dt_min_code).
- **Caveats**: MR is the NET flux through R_acc and includes the envelope's breathing (in the
  no-stream viper runs it swings by ~10 rho_stream Rsun^3 within 0.2 orbit, j/M then sits
  at the star's rotation 0.212); the stream mass over half an orbit is comparable, so read
  the j of accreted gas against that; cumulative integrals restart at zero on a restart
  (the script re-joins them); half an orbit is a feasibility test, not a steady state.

## 3. Back to viper (git, accretor-1006 only, never force)

```bash
cd <ak> || exit 1
git add docs/handover/plaskett-1007/results/numbers_*.json
echo "- $(date -u +%F) DeltaAI: Plaskett page <URL>; arms q12_s1, q12f_s1 t_final <x> orbits" \
  >> docs/handover/NOTE-2026-10-07-deltaai-plaskett.md
git add docs/handover/NOTE-2026-10-07-deltaai-plaskett.md
git commit -m "plaskett-1007: page URL + numbers JSONs"
git push fork accretor-1006   # if rejected: git fetch fork && git merge --no-edit fork/accretor-1006, push again
```

Fallback only if the Artifact tool is not available: also `git add` the PNGs in
`docs/handover/plaskett-1007/results/` (and the two `ryper.user.hst` files if < 5 MB each,
as `results/<arm>.user.hst`) and say so in the NOTE line.
