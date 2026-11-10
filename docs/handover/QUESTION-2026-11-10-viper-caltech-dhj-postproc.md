# QUESTION 2026-11-10 viper -> Caltech (USER): dhj round-1 post-processing on the viper/Raven runs: what we are missing

## Context

We are running your round-1 page pipeline on our own round-1 runs. Your tools are copied verbatim from
`docs/handover/caltech-to-raven-1107/pages/w121_dhj_r1/` into `/viper/ptmp2/jinma/dhjpost_1110/`. We changed only the
machine-specific bits:
- paths;
- RUNS;
- EOSF;
- sys.path;
- no r1b grey overlay;
- budget.SPECIAL emptied, because its text is about your TOI-1518b binary 2daa0a7a segment.

Our runs:
- WASP-76b, TOI-1518b, WASP-178b, WASP-189b and TOI-2109b: `dhj_r1_1110`, hydro 1x, prod-1110;
- WASP-121b hydro: `w121fresh_1107`, to rot 259;
- WASP-121b MHD 30 G dipole: continued from that hydro run at rot 259;
- WASP-121b FS0 reference.

The reductions, budgets, figures and cross-machine table (against your `dhj_p1_page_cache`) are done.

The items below exist only on Caltech. Following the user's rule, we do not rebuild them. Please push them to
`bsg-files-1009`, for example under `docs/handover/caltech-to-raven-1107/pages/w121_dhj_r1/`, or tell us how to run them.

## 1. Fe I fe_hires: the fast variant and the page hook

`fe_step.sh` calls the following, and none of them are on viper:
- `$J/fewind_1102/fehires_fast/fe_update.sh` (`--max 8 --img ...`);
- `fe_section.py`, which writes the page block and `fe_vs_rot.png`;
- `venv_fehgpu`.

Please push:
- `fehires_fast/`: fe_update.sh, fe_section.py and the GPU fe3d variant;
- the venv requirements (pip freeze of venv_fehgpu);
- which dumps `--max 8` selects (the last 8 dumps of each run?);
- whether `fe_vs_rot.png` averages over dumps.

Also confirm that fehires_fast gives the same number as the portable CPU tool
(`docs/handover/scripts/fe_hires`, cdb74ff28 Kp-marg) on the same dump. Your page quotes "Kp-marg (a)". We use the
`Kp-marg full (a) median [16, 84]` line of ana_summary.txt. Is that the same estimator?

## 2. The fe_hires planet JSONs you used for WASP-76b, TOI-1518b, WASP-178b and WASP-189b

The portable tool only ships `wasp121b*.json`. For the four planets we need the exact files you ran:
- R*, Rp, a/R*, b, P, Kp;
- LSF (the mixed ESPRESSO + MAROON-X choice for W76 and W189);
- nph, phase_range, rm_mask_phase;
- kp_marg_range;
- comptab 1x.

Until we have them, we ran only WASP-121b (hydro and MHD, with `wasp121b_1x_caltech.json`). TOI-2109b is not in
Seidel+2026. If you have a JSON for it, push it too.

## 3. obscheck_1102: phase curves, emission spectra, time-resolved Fe I

Your page section "Other observations" is static output of `$J/obscheck_1102` (em3d GPU RT, tr_fe). The code is not
on viper. Please push the scripts and their inputs (or their md5 and provenance if they are large), plus the command
lines that made `obs_*.png`. Without them we leave that section empty.

## 4. EOS table check

Our EOSF is a 1x AthenaK general-EOS dump:
- grid 381 x 451, log rho -19..0 step 0.05, log T 1.5..6.0 step 0.01;
- X 0.7381, Y 0.2485, H2 + ion + metal ion;
- md5 8d009459a9c8e68bcab854cb1f4c9aa6.

On its range (log rho >= -14) it is identical to viper's PLANET_EOSD 1x table, deepconv_0925.

Your dhjcs.py comment says your `anchor/eos_1x.txt` has 1126 T rows. Please send its header lines 1-6 and its md5,
so we can confirm both pages invert (rho, e) the same way.

Reply in a NOTE. The viper page is built without these items; each one is marked "awaits Caltech" on it.
