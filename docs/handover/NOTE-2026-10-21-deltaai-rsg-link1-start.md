# NOTE 2026-10-21 deltaai -> raven: RSG E link 1 STARTED on 8 GH200 (follow-up to NOTE-2026-10-21-deltaai-rsg-go-submitted)

- Link 1 = ghx4-interactive 3425380 on gh[072,137], started 2026-10-20 18:56 CDT (FRESH, t=0), 2 h link (-t 01:51:52).
  Log prints bin md5 c56693cb..., deck 9bd90b98..., template 4c442452... (all as reported).
- Setup ~150 s (cycle 0 at elapsed 151 s). Rate cycles 0->50: 2.25 s/cycle. dt = 1.0e3 s (flat so far, not the ~250 s in the
  deck comment). Latest: elapsed=2.635399e+02 cycle=50 time=5.000000e+04 dt=1.000000e+03
- GPU memory (nvidia-smi, per GH200, 8 ranks): ~52.5-52.9 GB of 96 GB (memcut OK; Raven measured 49 GB).
- Clean: 0 FATAL / DIVERGED / NaN. fofc fires only in x1 bin [546-547] (outer 2 cells), at-floor=0 ("other" ~6-8e5 cell-stages).
- At this rate a 2 h link is ~3000 cycles (~35 d of sim time if dt stays 1e3 s). rst every 0.5 d.
- SU: link 1 will be ~32 SU (interactive x2). Queue: ghx4 3425378 now estimated 2026-10-21 ~19:04 CDT (was 10-25), 3425379 afterok.
