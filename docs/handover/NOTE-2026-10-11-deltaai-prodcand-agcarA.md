# NOTE DeltaAI -> raven: AG Car A production-candidate gate (TASK-2026-10-11-deltaai-prodcand-agcarA) -- CAND DONE, ref pending

Picked up 2026-10-10 23:47 CDT. **cand1 + cand2 results posted 00:15 CDT; ref is re-queued (see below) and follows in this NOTE.**
Job **3357631** (ghx4-interactive gh113, 1 node x 4 GH200, 1 rank/GPU, MPICH_GPU_SUPPORT_ENABLED=1), wall 1 min 34 s.
- CAND prodcand-1011 e66ee669: `athena_hes_gpu_e66ee6698a6d` md5 c640dce072e336afa14da82054efa016.
- REF rt-integration d28dacd8: `athena_hes_gpu_d28dacd8aaa6` md5 315c33f73ffeffeeb7eb027c675d5e53.
- Per-run deck copies: nlim 60, ndiag 1, implicit_cfl 1251.3, output2/3/6/7 dcycle 60, rst block removed (instead of dt 1e30);
  cand: the 7 TASK keys, nothing else (face, fs_relax, h2be, hr_face sf, f0_init/rad_inertia/nofb_work = binary defaults).
- **ref FAILED at startup in 3357631** (rc 134): `Variable 'm1_thsw' in block 'output6' ... is not a valid choice`. REF d28dacd8 has
  no m1_thsw output; my mistake, I did not check REF against the deck's output variables. Re-queued as **3357651** with <output6>
  removed from the ref deck only (m1_fs kept: E/J_fs is available for ref); PENDING, priority 1820, estimated start 00:41 CDT.

## cand1 / cand2
| run | rc | FATAL | NON-CONV | Picard mean / max (solves) | s/cycle steady / rotation / average | cost, s per 1e4 s simulated (steady / average) | ic_rad start |
|---|---|---|---|---|---|---|---|
| cand1 | 0 | 0 | 0 | 3.25 / 7 (132) | 0.226 / 0.223 / 0.248 | 1.82 / 2.00 | 2.98 s (fs 0.32 s, relax 13 steps 2.23 s) |
| cand2 | 0 | 0 | 0 | 3.25 / 7 (132) | 0.227 / 0.224 / 0.245 | 1.83 / 1.97 | 2.96 s |

- **No rotation spike**: rotation cycles cost the same as steady ones (0.223 vs 0.226 s). (At the same 4-GH200 size, the earlier
  agcarA_1011 runs on vetacc-era code had rotation cycles of 0.45-1.0 s against ~0.20 s steady.)
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

## cand1 vs ref, ref timing, GPU-h
Pending on 3357651.

## Files
On DeltaAI: /work/nvme/bivj/jma20/prodcand_1011 (job.sub, analyse.py, runs/3357631/, logs/).
