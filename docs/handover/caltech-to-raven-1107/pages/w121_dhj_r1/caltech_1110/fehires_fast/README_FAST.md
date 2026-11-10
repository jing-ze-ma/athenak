# fehires_fast: fe_hires with one trace for both pressure cuts, and a GPU path (11-02)

`code/` is a copy of `../fehires/code`. Only `fe3d.py` was changed, and `fe3d_gpu.py` is new.
`ana.py`, `fe_core.py`, `csio.py` and `planets/` are unchanged. **The physics is unchanged.**

## What changed

- **`pcuts=1e-6,0`**: one ray trace makes both cuts.
  - The work that does not depend on the cut is done once: dump load, rays, KD-tree, interpolation, continuum, line profiles, exp(-tau_c).
  - The output is one `spec_<tag>_pc6.npz` and one `spec_<tag>_pc0.npz`. These are the same files and tags as the old one-cut runs, so `ana.py` and `collect.py` work unchanged.
  - Any short list of cuts works (up to 8), e.g. `pcuts=1e-6,1e-5,0` gives `_pc6`, `_pc5`, `_pc0`. Explicit names: `tags=a,b`.
  - `pcut=` alone (one cut) behaves as before.
- **`gpu=1`** (or env `FEHIRES_GPU=1`): the rays run on one GPU (CuPy, `code/fe3d_gpu.py`).
  - CPU worker processes (`NPROC`) do the ray geometry and the KD-tree lookup.
  - The GPU does the interpolation, the line-profile DGEMM, the line pruning, and the deposit + absorption.
  - The deposit + absorption is a per-pixel gather that keeps the CPU summation order, with no FMA.
  - Without `gpu=1` the CPU Pool path runs as before (fallback).
- **Venv:** `venv_fehgpu` (`/resnick/groups/carnegie_poc/jingze/fewind_1102/venv_fehgpu`).
  - Same numpy/scipy/numba/h5py as `venv_feh`, plus cupy-cuda12x 13.4.1.
  - No module load is needed.
  - `venv_feh` also runs the CPU path but has no CuPy.

## Run one dump on the GPU, both cuts

Take a free GPU first: `mkdir` lock in `/resnick/groups/carnegie_poc/jingze/lhllc2_1102/dhj`, and check that nvidia-smi shows < 100 MiB in use.

```
srun --overlap --jobid=<alloc> -N1 -n1 -c8 bash -c 'export CUDA_VISIBLE_DEVICES=<g>; \
  /resnick/groups/carnegie_poc/jingze/fewind_1102/fehires_fast/run_gpu.sh <planet> <dump.bin> <tag> <outdir> [eos_table]'
```

- `<planet>` is the json name in `code/planets` (w121, w76, w178, w189, toi1518).
- `eos_table` defaults to `w121r1_1101/anchor/eos_1x.txt`, the table that the 1x batch uses.
- Outputs, all in `<outdir>`:
  - `spec_<tag>_pc6.npz`, `spec_<tag>_pc0.npz`
  - `ana_<tag>_pc*.npz`, `ana_summary.txt`
  - logs
  - When `planets/<planet>_mx.json` exists, the MX ana goes into `<outdir>_mx`, as `run1.sh` does.
- The same thing by hand:
  ```
  cd fehires_fast/code; export OMP_NUM_THREADS=1 NPROC=8
  ../../venv_fehgpu/bin/python fe3d.py --planet planets/w121.json --eos EOS --out OUT DUMP tag=w121_00600 pcuts=1e-6,0 gpu=1
  ../../venv_fehgpu/bin/python ana.py --planet planets/w121.json --out OUT w121_00600_pc6   # and _pc0
  ```
- The usual fe3d keys still work (`phases=`, `nph=`, `dth=`, ...). For time-resolved runs, use `phases=a,b,c`.
- Cost (v3, 11-02 night): **about 100 s per dump for both cuts**, `-c 8` + 1 H200; ana adds about 10 s per cut.
  - Dump load: 2 s (`csio.EOS.invert` is now a numba scan, bitwise identical; it was 14 s).
  - GPU: about 3 s per exposure, 31 exposures. That splits into the line deposit + absorption gather 2.0 s, interpolation 0.34, profiles 0.2, DGEMM 0.34 and line pruning 0.13.
  - CPU prep (rays, KD-tree) runs in 8 workers, hidden behind the GPU.
  - GPU memory is a few GB (env `FEHIRES_GPU_BATCH`, default 256 rays per batch).
- The gather is the floor of the exact algorithm. Every ray deposits about 5,600 kept lines x 302 pixels and absorbs on all 731,800 pixels, for 36k rays x 2 cuts per exposure.
  - The CCF template is positive on 100 % of the pixels, so a template-pixel mask is not exact and saves nothing.

## Page series (fe_update.sh)

`fe_update.sh [--stride 10] [--max 8] [--img DIR] [--wait SEC]` keeps a per-planet blueshift-vs-rotation series in `series/`.

- It is incremental: the newest 5 dumps per planet, plus one dump every `stride` rotations. W121 and W76 stop at rot 300.
- It reuses `fehires/out`.
- It runs on a free test GPU (mkdir lock, owner fehgpu). If none is free, it falls back to one dump on the CPU one-pass path.
- It writes `series/img/fe_vs_rot.png`, `series/img/fehires_vs_seidel.png` (newest 5-dump window) and `series/img/fe_takeaway.txt`, and prints the paths.
- `python fe_section.py` prints the page's HTML section from `series/fe_series.json`.
- Spectra are deleted after ana; only the `ana_*.npz` files are kept.

## Validation (11-02; dumps W121 00600 and W178 00177, both cuts, against `../fehires/out`; v2 and v3 GPU and CPU one-pass)

- Spectra: **bitwise identical** float32 `spec` (0 of 45,375,444 pixels differ), for both the GPU and the CPU one-pass path.
- `ptan`: identical.
- `W` (contribution weights): agree to 7e-13 relative. Only the summation order differs.
- `ana`: Kp-Vsys map identical (max |diff| 0). Kp-marg, Kp0 peak and map maximum are identical.

| dump | cut | Kp-marg median full / rmcut (km/s) | map max (Kp, v) |
|---|---|---|---|
| W121 00600 | pc6 | -5.5939 / -5.5743 | 205.2, -4.0 |
| W121 00600 | pc0 | -5.8304 / -5.7509 | 205.2, -4.1 |
| W178 00177 | pc6 | -8.7877 / -8.1433 | 174.4, -8.9 |
| W178 00177 | pc0 | -8.9232 / -8.2286 | 174.4, -9.0 |

Validation outputs are in `val_gpu/`, `val_cpu/` and `val_old/`. The scripts are in `test/` (`cmp.py` does the comparison).
