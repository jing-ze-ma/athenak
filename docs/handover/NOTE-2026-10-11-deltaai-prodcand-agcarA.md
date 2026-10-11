# NOTE DeltaAI -> raven: AG Car A production-candidate gate (TASK-2026-10-11-deltaai-prodcand-agcarA) -- DONE

Picked up 2026-10-10 23:47 CDT. cand1 + cand2 posted 00:15 CDT; **ref added 00:30 CDT (job 3357651, gh147)**. All three runs rc 0, FATAL 0, NON-CONVERGED 0.
Job **3357631** (ghx4-interactive gh113, 1 node x 4 GH200, 1 rank/GPU, MPICH_GPU_SUPPORT_ENABLED=1), wall 1 min 34 s.
- CAND prodcand-1011 e66ee669: `athena_hes_gpu_e66ee6698a6d` md5 c640dce072e336afa14da82054efa016.
- REF rt-integration d28dacd8: `athena_hes_gpu_d28dacd8aaa6` md5 315c33f73ffeffeeb7eb027c675d5e53.
- Per-run deck copies: nlim 60, ndiag 1, implicit_cfl 1251.3, output2/3/6/7 dcycle 60, rst block removed (instead of dt 1e30);
  cand: the 7 TASK keys, nothing else (face, fs_relax, h2be, hr_face sf, f0_init/rad_inertia/nofb_work = binary defaults).
- **ref FAILED at startup in 3357631** (rc 134): `Variable 'm1_thsw' in block 'output6' ... is not a valid choice`. REF d28dacd8 has
  no m1_thsw output; my mistake, I did not check REF against the deck's output variables. Re-queued as **3357651** with <output6>
  removed from the ref deck only (m1_fs kept: E/J_fs is available for ref). 3357651 ran 00:23:49, wall 39 s, rc 0.

## Every run
| run | rc | FATAL | NON-CONV | Picard mean / max (solves) | s/cycle steady / rotation / average | cost, s per 1e4 s simulated (steady / average) | ic_rad start |
|---|---|---|---|---|---|---|---|
| cand1 | 0 | 0 | 0 | 3.25 / 7 (132) | 0.226 / 0.223 / 0.248 | 1.82 / 2.00 | 2.98 s (fs 0.32 s, relax 13 steps 2.23 s) |
| cand2 | 0 | 0 | 0 | 3.25 / 7 (132) | 0.227 / 0.224 / 0.245 | 1.83 / 1.97 | 2.96 s |
| ref | 0 | 0 | 0 | 3.78 / 9 (60) | 0.170 / 0.169 / 0.232 | 1.37 / 1.87 | - (default IC) |

- **Cost: cand is +33 % per steady cycle vs ref** (0.226 vs 0.170 s; hesdirk2 = 2 implicit stages, 132 solves vs 60, at
  lower Picard 3.25 vs 3.78) and **+7 % on the average** (0.248 vs 0.232 s/cycle, which includes the final-dump cycle).
- Rotation cycles: no spike in **either** binary here (cand 0.223, ref 0.169 s, = steady). **Correction to my 00:15 version**,
  which credited vgdrot for this: ref (no vgdrot) shows no spike either, so this test does not resolve the vgdrot gain at 4 GH200
  (Raven's BSG/AG Car A gpudev timings are the measurement). The 0.46-1.0 s rotation cycles I quoted came from icgen_1011 runs.
- ic_rad: `done, 2.983 s (fs 0.316 s, relax 13 steps 2.229 s), L_top = 2.8448341448e+38 (L/L_IC = 1.000612)`.
- **cand2 is bitwise vs cand1**: 12/12 bin files array-identical (hydro_w, m1, m1_thsw, m1_fs at cycles 0 and 60), user.hst and hydro.hst identical.

## L_top, cand1 (erg/s): flat from cycle 0
Range over cycles 0-60: 0.999689-1.000000 x L_top(0) (min at cycle 41). Every cycle:
0 2.844834e38, 1 2.844648, 2 2.844559, 3 2.844500, 4 2.844424, 5 2.844361, 6 2.844322, 7 2.844291, 8 2.844247, 9 2.844211,
10 2.844184, 11 2.844069, 12 2.844100, 13 2.844131, 14 2.844120, 15 2.844101, 16 2.844096, 17 2.844086, 18 2.844079, 19 2.844073,
20 2.844068, 21 2.844056, 22 2.844048, 23 2.844057, 24 2.844061, 25 2.844053, 26 2.844057, 27 2.844045, 28 2.844033, 29 2.844033,
30 2.844034, 31 2.844126, 32 2.844085, 33 2.844032, 34 2.844021, 35 2.844036, 36 2.844028, 37 2.844031, 38 2.844031, 39 2.844032,
40 2.844032, 41 2.843949, 42 2.843995, 43 2.844034, 44 2.844038, 45 2.844032, 46 2.844037, 47 2.844035, 48 2.844036, 49 2.844037,
50 2.844037, 51 2.844027, 52 2.844021, 53 2.844035, 54 2.844044, 55 2.844038, 56 2.844041, 57 2.844042, 58 2.844044, 59 2.844048, 60 2.844059 (x1e38).
Small steps at the rotation cycles (11, 31, 41: <= 3e-5 relative).

## E/J in the envelope at cycle 60, cand1 (median, p5-p95 over cells; cand2 identical)
| | 1.0-1.05 R_ph | 1.05-1.5 | 1.5-2 | 2-3 |
|---|---|---|---|---|
| E/J_fs (m1_e / m1_jfs) | 1.0024 (0.9862-1.0179) | **1.0008** (0.9715-1.0300) | 0.9991 (0.9235-1.1031) | 1.0064 (0.8397-1.2629) |
| E/J_tw (m1_e / m1_tj) | 1.0023 (1.0018-1.0028) | **1.0005** (0.9999-1.0015) | 1.0011 (1.0006-1.0045) | 1.0100 (1.0054-1.0157) |

ana_va.py is not in git, so these are my own: cellwise E/J over all cells of the shells in the band (R_ph 388.3 Rsun).
Same numbers as Raven's vetacc hpt arm at cycle 200 (1.0010 / 1.0008 in 1.05-1.5). The old production be/hr reading at
cycle 50 on DeltaAI (agcarA_1011 A runs) is E/J_fs 1.0194 / 1.0184 / 1.0076 / 1.0099.

## L_top, ref (for comparison)
ref: 2.843095e38 at cycle 0, **4.365807e38 at cycle 1 (x1.536)**, 4.107854 / 3.551379 / 3.150925 / 2.957766 at cycles 2-5,
2.847022e38 at 10, 2.844217e38 at 30, 2.844085e38 at 60 (x1.00035). Per-cycle table for both in deltaai-prodcand-1011/analysis.txt.
At cycle 60 cand1 / ref = 0.99999.

## E/J_fs in the envelope at cycle 60, ref (median, p5-p95)
1.0-1.05: 1.0197 (0.9976-1.0376) | 1.05-1.5: **1.0187** (0.9882-1.0474) | 1.5-2: 1.0078 (0.9316-1.1124) | 2-3: 1.0101 (0.8423-1.2658).
E/J_tw: n/a for ref (no thin switch, no m1_thsw output). So the 1.8 % E-over-J excess in 1.05-1.5 R_ph goes to 0.08 % (fs) / 0.05 % (tw) with cand.

## cand1 vs ref at cycle 60: max over shells in the band of |<q>_cand1 - <q>_ref| / |<q>_ref| (sin theta-weighted shell means)
| quantity | < 0.8 R_ph | 0.95-1.05 R_ph | > 1.05 R_ph |
|---|---|---|---|
| dens | 1.4e-5 | 2.95e-2 (1.004) | 1.54e-2 (1.193) |
| velx | 1.07 (0.406; velx ~ 0 there) | 9.2 (0.984; velx crosses 0) | 0.155 (2.799) |
| velx, max\|diff\| / max\|ref\| in band | 2.4e-2 (0.0003 of 0.014 km/s) | 4.9e-2 (0.28 of 5.70 km/s) | 3.2e-2 (0.29 of 8.84 km/s) |
| E | 1.2e-5 | 5.90e-2 (1.001) | 1.94e-2 (2.911) |
| F1 | 6.4e-3 (0.354) | 3.6e-5 | 1.1e-5 |

- The interior (< 0.8 R_ph) is unchanged to 1e-5 in dens and E; F1 differs by 0.6 % at the base (0.354 R_ph).
- The photosphere carries the change: E 6 %, dens 3 % at ~1.00 R_ph, consistent with ref still relaxing from its +54 % step-1 spike
  plus the new thin-envelope physics (E/J 1.019 -> 1.002).

## GPU-h
3357631 (cand1, cand2, failed ref): 1.6 min x 4 = 0.10; 3357651 (ref): 0.7 min x 4 = 0.04. **Total 0.15 GPU-h.**

## Files
In git: docs/handover/deltaai-prodcand-1011/ (job.sub, analyse.py, analysis.txt). On DeltaAI: /work/nvme/bivj/jma20/prodcand_1011 (runs/3357631, runs/3357651).
