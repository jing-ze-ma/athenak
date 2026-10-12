# NOTE DeltaAI: method-1012 fresh-run defaults bitwise check, full AG Car A 60 cycles (TASK-2026-10-12-deltaai-defaults) -- DONE: BITWISE (fresh defaults == all keys named on GH200)

**Picked up** 09:56 CDT. method-1012 **22f7021d** built: `athena_hes_gpu_22f7021dc1fb`, md5 9eaf7f298ab3d62c3ef28d407b11d6ee.
Decks resolved (@FILES@ -> agcarA_files copies, md5 4ab9e3f3 / b28e97c5 / 493f9d01), job.sub as provided (cmpout.py via venv python,
ATHENA_VIS = DeltaAI checkout vis/python); dry run OK (GU, then GN with @RV@ from GU's log).
Runs in ONE ghx4-interactive job **3366761** together with TASK rotslope-cost (defaults first, then rotslope FINC/FIN; 1 submitted job per user).
Dir /work/nvme/bivj/jma20/defaults_1012. STARTED line follows when it runs.

**STARTED**: job 3366761 RUNNING on gh143 since 10:18:55 CDT (defaults part first: GU, then GN, then cmpout).

## RESULT (posted 10:20 CDT): job 3366761, gh143, defaults part 10:19:01-10:20:30
**Verdict: fresh defaults == all keys named on GH200: yes.**
- cmpout.py: `GU vs .../runs/3366761/GN: 12 files compared, 0 differ -> BITWISE`
- GU log: `<rad_m1> vet_gd_rotate_dt = auto: 20 x the first dt 1252.56 = 25051.19312194301 (recorded; t0 0, align on, stretch 0.02)`
  -> GN deck line 99 `vet_gd_rotate_dt = 25051.19312194301`

| run | rc | FATAL | NON-CONV | Picard (126 solves) | s/cycle steady (c>=2, c%20!=0) | rot cycles 20, 40 | avg c>=2 | init | wall |
|---|---|---|---|---|---|---|---|---|---|
| GU | 0 | 0 | 0 | mean 2.794, max 6 | 0.2123 (n=55) | 0.274, 0.290 | 0.2147 | 21.7 s | 43.4 s |
| GN | 0 | 0 | 0 | mean 2.794, max 6 | 0.2136 (n=55) | 0.275, 0.291 | 0.2160 | 21.0 s | 42.0 s |

(timing from rotslope-1012/tana.py out.txt 2 20; the last step, which writes outputs, excluded). Both end at t = 7.515358e4 s, cycle 60.
Binary athena_hes_gpu_22f7021dc1fb md5 9eaf7f298ab3d62c3ef28d407b11d6ee. GPU-h: ~0.1 (1.5 min x 4 GH200).
