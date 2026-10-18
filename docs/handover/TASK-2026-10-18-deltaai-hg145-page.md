# TASK Raven -> DeltaAI (2026-10-18): build + PUBLISH the He giant production page from YOUR hg145 run (CPU only, 0 GPU SU)

USER 10-18: "delta ai is more progressed on the run, so you should ask delta ai to publish the url with its snapshots".
Star: 2.754 Msun He star (MESA r15140 14p5M_widebinary_profile14; "14p5M" = the 14.5 Msun PROGENITOR, not the star mass).

## 1. What
Run the page bundle docs/handover/hg145-page-1018/ (bsg-files-1009) on your production copy (0 -> 28 d now; refresh after each
link) and UPDATE the existing artifact **https://claude.ai/artifact/VKtffoTwjWw34E4SA5aHqa** (the Raven draft). Do NOT create a
new artifact: `Artifact action=read url=<that>` first, then publish with `url=<that>`, file publish/index.html, and
`files` = {"movie3d/hg145.mp4": publish/movie3d/hg145.mp4, "movie3d/hg145_poster.png": publish/movie3d/hg145_poster.png}.
Title stays "He Giant 2.75M Production".

## 2. How (one command; CPU node, not a GPU partition; no GPU SU)
```
cd <your checkout>/docs/handover/hg145-page-1018
RUN=/work/nvme/bivj/jma20/hg145_prod_1018/run  BINDIR=/work/hdd/bivj/jma20/hg145_prod_1018/bin \
OUT=/work/hdd/bivj/jma20/hg145_prod_1018/page  ICF=<checkout>/docs/handover/hg145-prod-1018/ic_hegiant_14p5M_mk.txt \
VISPY=<athenak>/vis/python  PY=<python3 with numpy+matplotlib>  PVPY=<pvbatch, or a python3 with `pip install vtk` (EGL/OSMesa)> \
FF=<ffmpeg with libx264>  SITE=DeltaAI  LAYOUT="1 node x 4 GH200 (ghx4), 4 MeshBlocks 704x64x64" \
LOGGLOB="run*.log"  JOBNAME=<your link job name>  NO_SBATCH=1  ./upd_page.sh
```
- Deps: numpy, matplotlib; VTK >= 9 with off-screen rendering (ParaView pvbatch or the pip vtk wheel); ffmpeg/libx264;
  AthenaK vis/python/bin_convert.py. Reads only the hst (hegiant.user.hst), run logs and the hydro_w + m1 bins; READ ONLY on the run.
- Cost: ~4 s per dump to reduce, ~1 s per frame to render (112 dumps for 28 d: ~10 min). Use a CPU node / CPU partition
  (or an interactive CPU allocation), not the login node for the full 0-28 d pass, and not a GPU partition.
- Incremental: re-running only reduces/renders new dumps. The movie colour scales are FIXED for the whole movie (v_r: 99.5 %
  of |v_r| on the side faces; T_eff: 0.5-99.5 % over all dumps); if the scale moves by > 5 % as the run grows, all frames
  re-render automatically (the caption states the values).
- Your hst has duplicate rows around 27.5 d (the restart from 00055.rst): mk_figs.py keeps the LAST occurrence of each time,
  so nothing to do, but check the history panel is continuous there.
- If your bins and hst carry different file stems than `hegiant.*`, tell Raven.

## 3. Check by eye before publishing
Every figure (history, profiles, photosphere maps, the two slice figures) and 3-4 movie frames (early, ~14 d, ~27 d, last).
Look in particular at: L_top/L_in drifting 1.00 -> 0.95 by ~27 d (your flag), the 1.03-1.10 R_ph shell above the photosphere
(on Raven, rho grows x30 at 115-119 Rsun with v_r +-5 km/s by 6.75 d), the single-column R_ph bump at phi 20.4 deg, lat 4.9 deg
(grows from 6 d on Raven), and the L_top/L_in spike to 1.17 at 6.08 d on Raven (Picard 3 -> 4 after it). Report whether
your copy shows the same (it should be the same setup; Raven reached t 9.1 d at TASK time).

## 4. Report
NOTE-2026-10-18-deltaai-hg145-page.md on bsg-files-1009: the URL (must be the one above), t covered, frame count, mp4 size,
index.html size (keep < ~10 MB), and the by-eye findings of sect. 3. Refresh the page (same URL) after each link.
