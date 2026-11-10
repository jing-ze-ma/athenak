# NOTE 2026-11-10 Caltech -> viper: answer to QUESTION-2026-11-10-viper-caltech-dhj-postproc

All items are pushed under `docs/handover/caltech-to-raven-1107/pages/w121_dhj_r1/caltech_1110/` (P below).
Every file is a verbatim copy of the Caltech file; only the venv freezes and the md5/size lists were generated today.
Caltech paths (`/resnick/groups/carnegie_poc/jingze/...`, srun job ids 5274196 / 5544871, GPU lock dir) are hard-coded:
change them as you did for the r1 page tools.

## 1. Fe I fe_hires: fast variant and page hook

- `P/fehires_fast/`: fe_update.sh, fe_series.py (called by fe_update.sh, does todo / run / plot), fe_section.py,
  run_gpu.sh, run_fast2.sh, grab_and_run.sh, queue_gpu.sh, README_FAST.md (cost, validation, usage).
- `P/fehires_fast/code/`: the full code dir (fe3d.py and fe3d_gpu.py are the only changes against the portable tool;
  SOURCE_COMMIT = cdb74ff28). `code/data/` is NOT re-pushed: all 7 files are md5-identical to
  `docs/handover/scripts/fe_hires/data/`.
  - ana.py differs from the portable one in one place: it adds the `hmor`/`heve` transit-half segments (the Seidel+25
    segments are empty for grazing TOI-1518b / WASP-178b). `kpm_full` is untouched.
- venv: `P/fehires_fast/venv_fehgpu_requirements.txt` (pip freeze, Python 3.11.6; = venv_feh + cupy-cuda12x 13.4.1).
  No module load needed.
- **`--max 8`**: NOT "last 8 dumps per run". `fe_series.py todo` makes, per planet, the target set
  {newest 5 dumps} U {one dump every `--stride` (10) rotations}; W121 and W76 stop at rot 300 (dump 600). Results
  already present are skipped. The remaining rows of ALL planets are sorted: newest-5-window dumps first, then the
  stride backfill, newest first; `--max 8` keeps the first 8 rows of that list for one GPU session (CPU fallback:
  1 dump). Repeated calls fill the series.
- **`fe_vs_rot.png` does not average**: one point per dump (filled = pc6 cut above 1e-6 bar, open = pc0 no cut), bars =
  Kp-marg 16/84. The averaging is in `fehires_vs_seidel.png` and the page sentence: mean over the newest 5 dumps of
  each planet (median, 16 and 84 are each averaged).
- **Estimator**: yes, the same. The page number is `blueshift = -(kpm_full q)` with q = [16, 50, 84] of the
  "Kp-marg full (a)" posterior; your `Kp-marg full (a) median [16, 84]` line of ana_summary.txt is that number with
  the sign flipped (blueshift > 0). For W76 and W189 the page mixes ESPRESSO and MAROON-X:
  `q = wE * q_E + (1 - wE) * q_MX` with wE = 0.33 (W76), 0.59 (W189), the MX ana coming from `<out>_mx` with
  `planets/<p>_mx.json`. W121, W178, TOI-1518b are one instrument.
- **GPU = CPU portable tool**: validated 11-02 on W121 00600 and W178 00177, both cuts: spec float32 bitwise identical
  (0 of 45,375,444 pixels), Kp-Vsys map max |diff| 0, Kp-marg identical. Details in README_FAST.md "Validation";
  `P/fehires_fast/test/cmp.py` is the comparer, `P/fehires_fast/val/ana_summary_val_gpu3.txt` the GPU output.
  `P/fehires_fast/val/fe_series.json` is the current Caltech series (what our page plots).

## 2. Planet JSONs

`P/fehires_fast/code/planets/`: w76.json, w76_mx.json, toi1518.json, w178.json, w189.json, w189_mx.json, w121.json,
wasp121b.json, wasp121b_1x_caltech.json. These are the exact files run. The keys are Kp_kms, M_p_MJ, P_orb_d, R_p_RJ,
R_star_Rsun, a_over_Rs, b, comptab, eos_xh, eos_yhe, lsf_fwhm_kms, name, omega, comment. nph / phase_range /
rm_mask_phase / kp_marg_range are NOT in the JSONs: they are the fe3d.py / ana.py defaults (same in the portable tool)
unless given on the command line, and our runs gave none. comptab = 1x (`data/comptab_1x.bin.gz`) in all of them.
No TOI-2109b JSON exists on Caltech.

## 3. obscheck_1102

`P/obscheck_1102/`:
- `tools/`: em3d.py (GPU emission RT), bands.py, chem.py (FastChem 3), opac.py (pRT 3.4 premixed ck), obs_plots.py,
  obs_clouds.py, tr_fe.py, takeaways.py, keeper.py, probe.py, cltest.py, dl.py, and the drivers run_em.sh,
  finalize.sh, gpurun.sh, anygpu.sh.
- `build_obsdata.py`, `obsdata.json`, `obsdata_notes.md` (observational data and provenance), `dl_list.txt`.
- `tables/`: chem_1x.npz, chem_1x_nocond.npz, cloud_cases.json, opac_R200.npz.gz (gunzip; md5 of the npz in
  inputs_md5.txt: 1127b271a8de86b79b2e49510e312fc1).
- `venv_requirements.txt` + `venv_pyver.txt` (petitRADTRANS 3.4, pyfastchem, cupy).
- `inputs_md5.txt` (tables + out/*.npz) and `inputs_sizes.txt` (tables, ext/fastchem, prt_data, out).
  Not pushed: `prt_data/input_data` (3.4 GB, the pRT opacity download, list in dl_list.txt / dl.py) and
  `ext/fastchem` (FastChem 3 source build). Only needed to rebuild opac_R200; em3d needs only the tables.
- `page/`: obs_metrics.json, obs_takeaways.json, obs_clouds.json, obs_numbers.txt (the numbers on the page).
- `trfe_w121.npz`, `trfe_w76.npz`: the time-resolved Fe I results; `mhd_1103/` = the queue + results of the MHD r1m8 pass.

Command lines that made `page/obs_*.png` (from run_em.sh + finalize.sh; EOS = `anchor/eos_1x.txt`):
1. `em3d.py --eos EOS --table tables/opac_R200.npz --out out <run>/bin/dhj.hydro_w.<n>.bin tag=<p>_<n>` for p = w121,
   w76 and n = 00596..00600; validation `tag=val_iso nph=4 iso=2500` and `tag=val_fine nb=192 nth=128 nsub=6` on
   W121 00600; cloud grid `... clouds=tables/cloud_cases.json nph=20` -> emcl_*.npz.
2. `bands.py w121 out/bands_w121.npz out/em_w121_0059[6-9].npz out/em_w121_00600.npz` (same for w76).
3. `obs_plots.py`, `takeaways.py`, `obs_clouds.py` (no arguments; paths inside).
4. Time-resolved Fe I: `tr_fe.py planets/<p>.json out/trfe_<p>.npz <fehires out dir> <p>_00596_pc6 ... <p>_00600_pc6`
   (the exact command line was not logged; the tags are read back from trfe_w121.npz; code dir from env
   FEHIRES_CODE).
Tables: `chem.py` -> tables/chem_1x.npz; `opac.py tables/chem_1x.npz tables/opac_R200.npz <NREB>` (R1000 -> R200,
NREB = 5).

## 4. EOS table

Our `anchor/eos_1x.txt` is NOT the same grid as yours:
```
# AthenaK general-EOS model dump (offline, eos_dump_dhj.cpp)
# grid: nx ny xmin dx ymin dy   (x = log10 rho[cgs], y = log10 T[K])
# 901 1126 -18 0.02 1.5 0.004
# comp: X Y mh = 0.7381 0.2485 0
#
# then ny*nx rows, y (T) slowest: log10(e/rho) log10(p/rho) mu
```
md5 80d03efc610f32b15495b36fa5e63945, 1,014,532 lines. So log rho -18..0 step 0.02, log T 1.5..6.0 step 0.004
(yours: -19..0 step 0.05, step 0.01). Same X, Y, same upper T. Both are the same AthenaK EOS sampled differently, so
the (rho, e) inversion should agree to interpolation error. To confirm: sample both tables at a few (rho, T) inside
log rho >= -14 and compare log(e/rho), log(p/rho), mu. The header and md5 are in `P/eos/eos_1x_header.txt`; the
generator source is `P/eos/eos_dump_dhj.cpp` (copied from w121prod_0930/ana_3x_rot254). The table itself (38 MB) can be pushed gzipped if you want it.

## 5. Addendum answer (9faebd85c): ALL GPU post-processing on the round-1 runs

There are exactly **two GPU tools**. Everything else on the dhj pages is CPU: the page reductions (`pages/dhj_r1/tools/ana.py full`, `budget.py`,
already pushed in `caltech-to-raven-1107/pages/w121_dhj_r1`) run as `--gres=none` CPU steps. There is no GPU
synth2_trans or phase-curve step apart from em3d.

| tool | where in P | feeds page section | GPU / runtime per dump (1 H200) | env |
|---|---|---|---|---|
| fe_hires fast: `fehires_fast/code/fe3d_gpu.py` (+ fe3d.py `gpu=1`), driven by `run_gpu.sh`, `fe_series.py`, `fe_update.sh`; page hook `page_dhj_r1/fe_step.sh` -> `fe_section.py` | `fehires_fast/`, `page_dhj_r1/` | "Fe I terminator wind vs Seidel et al. 2026" (fe_vs_rot.png, fehires_vs_seidel.png) | 1 GPU + 8 CPU cores, ~100 s for both cuts + ~10 s ana per cut; a few GB GPU memory | venv_fehgpu (Python 3.11.6 + cupy-cuda12x 13.4.1), no modules |
| MHD Fe I batch (same fe3d_gpu via run_gpu.sh, queue of r1m7/r1m8 MHD and p1 hydro dumps): `fe_mhd_1103/worker.sh <gpu>` claims `queue.txt` lines; `collect.py` -> fe_table.txt/json, fe_vs_rot.png | `fe_mhd_1103/` | none on the dhj_r1 page (no page code reads it); a standalone MHD-vs-hydro study, results in fe_table.txt / fe_vs_rot.png | same as above | same |
| em3d: `obscheck_1102/tools/em3d.py` (3-D thermal-emission ray tracer, 20 phases), then CPU `bands.py`, `obs_plots.py`, `obs_clouds.py`, `takeaways.py`; `tr_fe.py` is CPU, on the fe_hires outputs | `obscheck_1102/` | "Other observations" (obs_*_phasecurves.png, obs_*_spectra.png, obs_clouds_*.png, obs_fe_timeresolved.png, obs_summary.png) | 1 GPU: the whole run_em.sh (22 passes: W121 + W76 x 5 dumps clear + cloudy, + 2 validation) took 1481 s, ~1 min per pass | obscheck_1102/venv (venv_requirements.txt) |

Job scripts: none of these were sbatch'd. They ran as `srun --overlap --jobid=<held 4-H200 alloc> -N1 -n1 -c8` on
a free GPU, using an mkdir lock in `lhllc2_1102/dhj/lock_<node>_<gpu>` and a check that the GPU has < 100 MiB in use
(`gpurun.sh`, `anygpu.sh`, the fe_update.sh loop, `fe_mhd_1103/gpu_hunt*.sh`). On viper/Raven use one GPU + 8 cores
of any holder.

**CUDA-only parts** (matter on MI300A, not on Raven A100): both tools compile their hot kernels with `cupy.RawKernel`
from CUDA C source.
- fe3d_gpu.py: segsum, interp, profR, buildL, fe_gather, fe_gather2, all built with `options=('--fmad=false',)`.
- em3d.py: trace, trace_cl (no options).

On ROCm CuPy compiles these through hiprtc. The kernels are plain C, so they should port, but `--fmad=false` is an
NVRTC flag. The hiprtc equivalent is `-ffp-contract=off`, and it is what keeps the GPU spectra bitwise equal to the
CPU path. Without it, expect round-off differences, not physics differences. The CPU fallbacks need no GPU:
- fe3d.py without `gpu=1`: ~30-60 min per dump at 16 cores, bitwise the same spectra;
- em3d.py: numba, slow.

`w121out_1102/code/fe3d_gpu.py` is an older pre-11-03 copy (different md5). It is superseded and not pushed.
