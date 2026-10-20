# TASK 2026-10-20 raven -> caltech: build + PUBLISH the 3-D accretor page from YOUR fresh tracer run (CPU only, 0 GPU SU)

USER 10-20: Caltech updates the 3-D accretor page with its fresh tracer run (acc3d_trc_1019, deck v11amwtr, binary prodmerge-1019b
13206541, 2 H200). Raven's run11 (stopped at t 0.73) stays only as a "previous run" reference (fixed images in the bundle).

## 1. What
Run the portable bundle `docs/handover/acc3d-page-1020/` (bsg-files-1009) on your run, READ ONLY on the run, and UPDATE the
existing artifact **https://claude.ai/artifact/U37ZLacbA1JBKRaFz2MzQa** ("RY Per Accretor 3-D"). Do NOT create a new artifact.

## 2. Dependencies
- python3 >= 3.9 with numpy + matplotlib (your page modules work: `module load python/3.11.6-gcc-11.3.1-rphh4kv py-matplotlib/3.8.2-gcc-11.3.1-7d3v7ui`).
  No h5py / scipy / AthenaK checkout needed: the bundle carries `binread.py` (AthenaK `read_binary`, format 1.2 with the true
  stretched faces mb_x1f/mb_x2f/mb_x3f). VISPY=<athenak>/vis/python uses that bin_convert.py instead (it must have read_binary with mb_x1f).
- ffmpeg with libx264 (your spack ffmpeg-6.0 path from bsgprod-1012/pagework/make_page.sh works).
- VTK >= 9 with off-screen rendering for the 3-D movie: ParaView `pvbatch`, or a venv python with `pip install vtk` (EGL) or
  `pip install --extra-index-url https://wheels.vtk.org vtk-osmesa` (no GPU/X needed). render3d.py is pure VTK (no paraview.simple).
  If you have neither: SKIP3D=1 (the page then shows only the 2-D movie; say so in the NOTE).
- The run's EOS table dump (25.7 MB, too big for git): `eos_x0.70_y0.28_logd-18.txt`, md5 **e548d9c2777b8b63f979c024ab2dd564**
  (Raven: /raven/ptmp/jinma/acc3d_1017/mk/eosdump/). Either ask Raven to scp it, or regenerate it with your binary
  (it is written at startup; 0 cycles): `BN=<your 13206541 binary> DK=<your substituted deck acc3d_prod_trc.athinput> DEST=<dir>
  LAUNCH="srun -n 1" ./mk_eostbl.sh` inside a 1-GPU allocation (~2 min), then check the md5. Pass it as EOSTBL. The bundle warns
  (does not stop) on a different md5; report it.
- The Rosseland table of the deck for the optical depths: default `../acc3d-prod-1018/rosseland_mesakap_gs98_x0.70_z0.02.txt`
  in the same checkout (md5 d81b7a46, = your files/ copy). Override with OPAC.

## 3. How (one command; CPU node, e.g. `sbatch -A carnegie_poc -p expansion -c 8 --mem=64G -t 04:00:00`)
```
cd <your bsg-files-1009 checkout>/docs/handover/acc3d-page-1020
RUN=/resnick/groups/carnegie_poc/jingze/acc3d_trc_1019 \
OUT=/resnick/groups/carnegie_poc/jingze/acc3d_trc_1019_page \
EOSTBL=/resnick/groups/carnegie_poc/jingze/acc3d_trc_1019/files/eos_x0.70_y0.28_logd-18.txt \
PY=python3 FF=<ffmpeg with libx264> PVPY=<pvbatch | venv python with vtk> \
SITE=Caltech LAYOUT="1 node x 2 H200 (hpc-sm), 8 MeshBlocks 420x58x64 (4 per GPU); 7-d links" \
BIN_LABEL="prodmerge-1019b @ 13206541 (H200 nofma md5 201918b0)" DECK_LABEL="acc3d_prod_v11amwtr (substituted md5 9d1b56bb)" \
RUN_LABEL="fresh tracer run acc3d_trc_1019" CHAIN_LABEL="acc3d_chain.sh, 7-d links 5088058-63, fresh from t = 0" \
JOBNAME=<your link job name, for the queue line> NP=8 ./upd_page.sh
```
- RUN is the parent of static/ and stream/ (the bundle reads ryper.{user,hydro,trc}.hst, run.log, bin/ryper.{hydro_w,m1}.NNNNN.bin,
  ../chain.txt; nothing is written there). Everything else goes to OUT (data caches, img/, movie frames, publish/).
- Off Raven the script never calls sbatch; it runs in place (put it in your own batch job as above).
- Incremental: caches per dump in OUT/data; 2-D frames only for new hydro dumps that have an m1 dump at the same time (every 0.01),
  3-D frames for every new hydro dump, tracer analyses for NTRC (4) times spread over the stream phase + the latest dump.
  FORCE=1 redraws everything.
- Expected runtime (Raven CPU, 8 cores, NP 6, prod-size dumps): see sect. 6; roughly 25 s per new hydro_w dump for the 3-D
  reduction + frame, ~20 s per paired dump for the 2-D frame, ~1 min per tracer-analysis dump, plus ~1 min fixed (EOS table parse,
  cached after the first run). The first full pass on ~150 dumps is ~1-1.5 h at NP 8; later refreshes only process new dumps.
  Memory ~2 GB per worker.

## 4. Check by eye before publishing
- Header: t, phase, site/layout, binary, deck must describe YOUR run (env labels above).
- Run checks: NaN/FATAL/NON-CONVERGED all 0; dt retries; the mass budget residual ~1e-9 relative; energy residual vs E_bnd
  (accwall gate 1e-3); the θ-wall screen.
- Tracers: dNz0/dt vs 254.84 (you saw 254.9); Ms0 budget residual at round-off; s0 maps / profiles; the s0-weighted entropy table
  (run11's top-dM estimate was s_gas 38.5, p10-90 38.0-39.1; the tracer value replaces it); the Ω(r) plot.
- Maps + 2-D movie: columns log ρ, log T_gas (EOS table), log P_rad/P_gas (P_rad = E/3 from m1, P_gas from the EOS table), v_φ.
  The thin ambient is radiation dominated (P_rad/P_gas up to 1e5, saturates the red end at 1e5); envelope ~0.1-10.
- 3-4 movie frames (early static, contact ~t 0.068, mid, last).

## 5. Publish (same URL)
`Artifact action=read url=https://claude.ai/artifact/U37ZLacbA1JBKRaFz2MzQa` first, then publish with that `url`,
`file_path` = OUT/publish/index.html and
`files` = {"movie/acc3d_eq.mp4": OUT/publish/movie/acc3d_eq.mp4, "movie/acc3d_eq_poster.png": OUT/publish/movie/acc3d_eq_poster.png,
"movie/acc3d_3d.mp4": OUT/publish/movie/acc3d_3d.mp4, "movie/acc3d_3d_poster.png": OUT/publish/movie/acc3d_3d_poster.png}
(the page references exactly these relative paths). Title stays "RY Per Accretor 3-D". index.html is ~4-6 MB (figures inlined);
keep it < ~10 MB. Refresh after every link (and once a day while a 7-d link runs), same URL.

## 6. Raven test of the bundle (10-20)
See the NOTE for this TASK's commit: bundle run end to end on Raven with env changes only, on (a) a run11 subset (11 hydro +
8 m1 dumps; no trc.hst, no s_00: tracer section skipped cleanly) and (b) the acctrc gate g19c2 C + Sm (trc.hst + s_00/s_01).

## 7. Report
NOTE-2026-10-20-caltech-acc3d-page.md on bsg-files-1009: URL (must be the one above), t covered, dump / frame counts, mp4 sizes,
index.html size, runtime, EOS table md5, and the by-eye findings of sect. 4 (in particular the s0-weighted entropy and where the s0 gas sits).

## 8. Known issues (ACCPAGE review 10-20; NOT fixed, USER: push as is)
- mk_eostbl.sh runs the 0-cycle EOS dump with `meshblock/nx3=512`, i.e. ONE 420x58x512 MeshBlock (12.5M cells) on one GPU, which
  holds 2x the per-GPU cells of the production layout (4 x 420x58x64 per H200) and may run out of GPU memory. If it does, run it with
  `mesh/nx3=64 meshblock/nx3=64` instead (the table does not depend on the grid; untested), or simply scp the Raven copy
  (/raven/ptmp/jinma/acc3d_1017/mk/eosdump/eos_x0.70_y0.28_logd-18.txt, md5 e548d9c2777b8b63f979c024ab2dd564 verified 10-20).
- Map figures (mk_prod.py): the bottom-row titles "merid. through impact: ... (φ N)" overflow the panel in the P_rad/P_gas and v_φ
  columns (the "(φ N)" suffix is clipped). Cosmetic only.
- Tracer figures with no s0 gas on the star yet (early stream): the τ panel of fig_trc_entropy.png is empty with an empty legend box. Expected.
- The Raven end-to-end retest of the bundle was cancelled (USER 10-20: push now); the 10-20 Raven tests of sect. 6 are the earlier
  worker's (the run11-subset test log stops at the render3d step, the g19c2 trc test completed: index.html 1.5 MB, both mp4s built).
