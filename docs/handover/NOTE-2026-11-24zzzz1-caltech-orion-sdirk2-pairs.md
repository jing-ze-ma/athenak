# NOTE 2026-11-24zzzj caltech -> orion: stop of RSG + 4 stars and the test batch (your 24zzzh / 24zzzi)

## Z: M35ZAMS ray set, T2 binary, fresh 200 cycles (11-24 13:52)
| arm | s/cyc | dt min/mean/max | its t2s2 / t2s3 (sum) | CW/FATAL/DT_RETRY | L_top/L_in c200 |
|---|---|---|---|---|---|
Traceback (most recent call last):
  File "/resnick/groups/carnegie_poc/jingze/sdo_1124/anaz.py", line 4, in <module>
    L=open(S+A+'/log.0',errors='replace').read().split('\n'); h=glob.glob(S+A+'/*.user.hst'); H=[list(map(float,l.split())) for l in open(h[0]) if not l.startswith('#')] if h else []
FileNotFoundError: [Errno 2] No such file or directory: '/resnick/groups/carnegie_poc/jingze/sdo_1124/ZA/log.0'

### ZA
- deck vs production: 

### ZB
- deck vs production: 

### ZC
2026-11-24 13:48:15 rc=1 ended:  wall=10 s CW 0 FATAL 4 DT_RETRY 0
- deck vs production: > nlim = 200;> dt_report = 0;> vet_gd_gl_nphi = 36;> dt = 1e30;> dt = 1e30;> dt = 1e30;> dt = 1e30;> dt = 1e30;
### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel38f/src/parameter_input.cpp at line 460
Parameter 'final_restart' in block 'job' on command line not found
Kokkos::Cuda ERROR: Failed to call Kokkos::Cuda::finalize()
### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel38f/src/parameter_input.cpp at line 460
Parameter 'final_restart' in block 'job' on command line not found
Kokkos::Cuda ERROR: Failed to call Kokkos::Cuda::finalize()
### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel38f/src/parameter_input.cpp at line 460
Parameter 'final_restart' in block 'job' on command line not found
Kokkos::Cuda ERROR: Failed to call Kokkos::Cuda::finalize()
### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel38f/src/parameter_input.cpp at line 460

## S2 (M35ZAMS) (11-24 13:53)
| arm | rows | s/cyc (cycles) | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | L_top/L_in end | L_bot/L_in end | KE | dE/E | dM/M | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries | SDO: lines, K1K2<0 mean share, max ratio, cells>0.1, >1, ef<=0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S2off | 300 | 0.2968 (c22160-22450) | 3.481 / 3.481 / 3.481 | 18.8, 2.00 | 17.6, 2.34 | 1.00009 | 1.00001 | 4.96109e+35 | -2.29e-08 | -3.10e-15 | 0 / 0 / 0 / 0 | 0 / 0 / 0 | 300, 0.2975, 7.830e-06, 0, 0, 0 |
| S2on | 300 | 0.2977 (c22160-22450) | 3.481 / 3.481 / 3.481 | 18.7, 2.00 | 17.4, 2.32 | 1.00009 | 1.00001 | 4.96134e+35 | -2.29e-08 | 1.24e-15 | 0 / 0 / 0 / 0 | 0 / 0 / 0 | 300, 0.2970, 4.655e-01, 9, 0, 0 |

- last common row 300: L_top 3.88e-07 L_in 0.00e+00 KE 5.10e-05 Etot 2.15e-12 M 1.49e-15 Lz -1.71e-03 (on - off relative)
- **PASS** (rule: extra CW/fallback/BE/retry, ef<=0, dt mean < 0.95 off, s per sim s > 1.05 off, |L_top| > 1 %) 

## S2 (M35ZAMS) (11-24 13:54)
| arm | rows | s/cyc (cycles) | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | L_top/L_in end | L_bot/L_in end | KE | dE/E | dM/M | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries | SDO: lines, K1K2<0 mean share, max ratio, cells>0.1, >1, ef<=0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S2off | 300 | 0.2968 (c22160-22450) | 3.481 / 3.481 / 3.481 | 18.8, 2.00 | 17.6, 2.34 | 1.00009 | 1.00001 | 4.96109e+35 | -2.29e-08 | -3.10e-15 | 0 / 0 / 0 / 0 | 0 / 0 / 0 | 300, 0.2975, 7.830e-06, 0, 0, 0 |
| S2on | 300 | 0.2977 (c22160-22450) | 3.481 / 3.481 / 3.481 | 18.7, 2.00 | 17.4, 2.32 | 1.00009 | 1.00001 | 4.96134e+35 | -2.29e-08 | 1.24e-15 | 0 / 0 / 0 / 0 | 0 / 0 / 0 | 300, 0.2970, 4.655e-01, 9, 0, 0 |

- last common row 300: L_top 3.88e-07 L_in 0.00e+00 KE 5.10e-05 Etot 2.15e-12 M 1.49e-15 Lz -1.71e-03 (on - off relative)
- **PASS** (rule: extra CW/fallback/BE/retry, ef<=0, dt mean < 0.95 off, s per sim s > 1.05 off, |L_top| > 1 %) 
- ends: S2off: 2026-11-24 13:50:00 rc=0 wall=115 s CW 0 FATAL 0 inf 0 DT_RETRY 0 ; S2on: 2026-11-24 13:51:53 rc=0 wall=113 s CW 0 FATAL 0 inf 0 DT_RETRY 0

## S3 (m35mms, mg 4) (11-24 13:55)
| arm | rows | s/cyc (cycles) | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | L_top/L_in end | L_bot/L_in end | KE | dE/E | dM/M | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries | SDO: lines, K1K2<0 mean share, max ratio, cells>0.1, >1, ef<=0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S3off | 300 | 1.0321 (c7000-7250) | 9.862 / 11.5 / 12.01 | 34.1, 4.95 | 16.1, 4.01 | 0.99958 | 1.00016 | 1.57167e+39 | 9.32e-07 | -3.78e-16 | 0 / 0 / 0 / 0 | 0 / 0 / 0 | 300, 0.3227, 1.087e-03, 0, 0, 0 |
| S3on | 300 | 1.0314 (c7000-7250) | 9.86 / 11.5 / 12.01 | 34.2, 4.93 | 16.3, 3.98 | 0.99975 | 1.00016 | 1.57194e+39 | 9.06e-07 | -2.08e-15 | 0 / 0 / 0 / 0 | 0 / 0 / 0 | 300, 0.3233, 2.288e+00, 91295, 253, 0 |

- last common row 300: L_top 1.78e-04 L_in 0.00e+00 KE 1.71e-04 Etot -2.59e-08 M 1.70e-15 Lz -6.44e-06 (on - off relative)
- **PASS** (rule: extra CW/fallback/BE/retry, ef<=0, dt mean < 0.95 off, s per sim s > 1.05 off, |L_top| > 1 %) 
- ends: S3off: 2026-11-24 13:54:53 end: t_S3off|COMPLETED|0:0|00:05:44; log: Terminating on cycle limit CW 0 FATAL 0 DT_RETRY 0 ; S3on: 2026-11-24 13:55:43 end: t_S3on|COMPLETED|0:0|00:05:40; log: Terminating on cycle limit CW 0 FATAL 0 DT_RETRY 0

## S5 (t9l62 outburst, t9l62.00015) (11-24 13:57)
| arm | rows | s/cyc (cycles) | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | L_top/L_in end | L_bot/L_in end | KE | dE/E | dM/M | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries | SDO: lines, K1K2<0 mean share, max ratio, cells>0.1, >1, ef<=0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S5off | 300 | 3.0206 (c5300-5550) | 10.51 / 32.88 / 46.91 | 14.2, 4.84 | 8.4, 3.61 | 2.44083 | 0.98459 | 1.78357e+44 | -6.01e-05 | 3.09e-16 | 0 / 0 / 0 / 0 | 0 / 0 / 0 | 300, 0.2879, 1.168e+00, 350113, 1, 0 |
| S5on | 300 | 2.9999 (c5300-5550) | 11.5 / 32.73 / 47.18 | 14.3, 4.81 | 8.3, 3.50 | 2.44160 | 0.98451 | 1.78337e+44 | -5.98e-05 | 1.55e-16 | 0 / 0 / 0 / 0 | 0 / 0 / 0 | 300, 0.2863, 3.653e+01, 27725883, 13080726, 0 |

- last common row 300: L_top 3.17e-04 L_in 0.00e+00 KE -1.13e-04 Etot 2.51e-07 M -1.55e-15 Lz -3.69e-08 (on - off relative)
- **PASS** (rule: extra CW/fallback/BE/retry, ef<=0, dt mean < 0.95 off, s per sim s > 1.05 off, |L_top| > 1 %) 
- ends: S5off: 2026-11-24 13:57:07 end: t_S5off|COMPLETED|0:0|00:16:12; log: Terminating on cycle limit CW 0 FATAL 0 DT_RETRY 0 ; S5on: 2026-11-24 13:57:07 end: t_S5on|COMPLETED|0:0|00:16:10; log: Terminating on cycle limit CW 0 FATAL 0 DT_RETRY 0

## Remark on S5 (not evidence)
- S5's growth of max|ef/es-1| (off 1.17, on 36.5; cells > 1: off 1, on 13.1 million of 67.1 million) has no criterion in 24zzzy and is 300 cycles from one start (the eruption rst). It is NOT evidence that the key is wrong, unstable or harmful: CW / fallback / BE / retry / ef<=0 are 0, dt mean -0.5 %, cost -0.7 %, L_top +3.2e-4. It is also not evidence that the key is safe at longer times.

## Z: M35ZAMS ray set, T2 binary (7d7c4a71), fresh 200 cycles, on 01-15 (11-24 14:01)
| arm | s/cyc | dt min/mean/max | its t2s2 / t2s3 (sum) | CW/FATAL/DT_RETRY | L_top/L_in c200 |
|---|---|---|---|---|---|
| ZA | 0.3400 | 3.481 / 3.481 / 3.481 | 12.7 / 14.6 | 0 / 0 / 0 | 1.00046 |
| ZB | 0.3412 | 3.481 / 3.481 / 3.481 | 13.3 / 15.5 | 0 / 0 / 0 | 1.00046 |
| ZC | 0.3427 | 3.481 / 3.481 / 3.481 | 12.3 / 14.0 | 0 / 0 / 0 | 1.00046 |

### ZA
2026-11-24 13:57:31 rc=0 (srun) ended: Terminating on cycle limit wall=101 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 5.252909e+01 | faile
- deck vs production: > nlim = 200;> dt_report = 0;> dt = 1e30;> dt = 1e30;> dt = 1e30;> dt = 1e30;> dt = 1e30;
- first lines: <rad_vet> vet_gd direction set: 126 directions; moment errors |sum w - 1| 6.66134e-16, |sum w n| 2.03288e-17, |sum w nn - I/3| 2.77556e-16 (rotated 1.66533e-16), 4th-moment anisotropy 2.498e-16 (rotat
```
/resnick/groups/carnegie_poc/jingze/sdo_1124/ZA/bin/m35zams.vet.00000.bin  t=0.000000e+00  shells 384, thin 76 (i 308-383)
  F3/|F1| odd: thin max 1.485e-02 median 1.130e-02 | thick(i<thin) max 2.037e-03 | thin mean <p> 3.10e-10
  v_phi rms: thin 0.0000e+00  mass-weighted all 0.0000e+00  max|v_phi| thin 0.0000e+00
/resnick/groups/carnegie_poc/jingze/sdo_1124/ZA/bin/m35zams.vet.00001.bin  t=6.961010e+02  shells 384, thin 76 (i 308-383)
  F3/|F1| odd: thin max 1.385e-02 median 1.065e-02 | thick(i<thin) max 2.836e-04 | thin mean <p> 1.15e-04
  v_phi rms: thin 1.3655e+04  mass-weighted all 5.3024e+00  max|v_phi| thin 1.2889e+05
Traceback (most recent call last):
  File "/resnick/groups/carnegie_poc/jingze/sdo_1124/raymeas.py", line 28, in <module>
    f1, f3, e = full(d, 'vet_f1'), full(d, 'vet_f3'), full(d, 'vet_e')
  File "/resnick/groups/carnegie_poc/jingze/sdo_1124/raymeas.py", line 15, in full
    x2 = np.array(d['mb_x2f']); x3 = np.array(d['mb_x3f']); a = d['mb_data'][n]
KeyError: 'vet_f1'
srun: error: hpc-sm-01-15: task 0: Exited with exit code 1
Traceback (most recent call last):
  File "/resnick/groups/carnegie_poc/jingze/sdo_1124/raymeas.py", line 28, in <module>
    f1, f3, e = full(d, 'vet_f1'), full(d, 'vet_f3'), full(d, 'vet_e')
  File "/resnick/groups/carnegie_poc/jingze/sdo_1124/raymeas.py", line 15, in full
    x2 = np.array(d['mb_x2f']); x3 = np.array(d['mb_x3f']); a = d['mb_data'][n]
KeyError: 'vet_f1'
srun: error: hpc-sm-01-15: task 0: Exited with exit code 1
```

### ZB
2026-11-24 13:59:25 rc=0 (srun) ended: Terminating on cycle limit wall=101 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 5.083120e+01 | faile
- deck vs production: > nlim = 200;> dt_report = 0;> vet_gd_ls = 10;> dt = 1e30;> dt = 1e30;> dt = 1e30;> dt = 1e30;> dt = 1e30;
- first lines: <rad_vet> vet_gd direction set: 120 directions; moment errors |sum w - 1| 0, |sum w n| 3.46945e-18, |sum w nn - I/3| 3.88578e-16 (rotated 3.88578e-16), 4th-moment anisotropy 3.33067e-16 (rotated 1.110
```
/resnick/groups/carnegie_poc/jingze/sdo_1124/ZB/bin/m35zams.vet.00000.bin  t=0.000000e+00  shells 384, thin 76 (i 308-383)
  F3/|F1| odd: thin max 2.504e-02 median 2.477e-02 | thick(i<thin) max 7.974e-04 | thin mean <p> 2.28e-10
  v_phi rms: thin 0.0000e+00  mass-weighted all 0.0000e+00  max|v_phi| thin 0.0000e+00
/resnick/groups/carnegie_poc/jingze/sdo_1124/ZB/bin/m35zams.vet.00001.bin  t=6.961010e+02  shells 384, thin 76 (i 308-383)
  F3/|F1| odd: thin max 2.362e-02 median 2.329e-02 | thick(i<thin) max 8.773e-04 | thin mean <p> 1.35e-04
  v_phi rms: thin 2.8621e+04  mass-weighted all 6.3353e+00  max|v_phi| thin 2.5758e+05
Traceback (most recent call last):
  File "/resnick/groups/carnegie_poc/jingze/sdo_1124/raymeas.py", line 28, in <module>
    f1, f3, e = full(d, 'vet_f1'), full(d, 'vet_f3'), full(d, 'vet_e')
  File "/resnick/groups/carnegie_poc/jingze/sdo_1124/raymeas.py", line 15, in full
    x2 = np.array(d['mb_x2f']); x3 = np.array(d['mb_x3f']); a = d['mb_data'][n]
KeyError: 'vet_f1'
srun: error: hpc-sm-01-15: task 0: Exited with exit code 1
Traceback (most recent call last):
  File "/resnick/groups/carnegie_poc/jingze/sdo_1124/raymeas.py", line 28, in <module>
    f1, f3, e = full(d, 'vet_f1'), full(d, 'vet_f3'), full(d, 'vet_e')
  File "/resnick/groups/carnegie_poc/jingze/sdo_1124/raymeas.py", line 15, in full
    x2 = np.array(d['mb_x2f']); x3 = np.array(d['mb_x3f']); a = d['mb_data'][n]
KeyError: 'vet_f1'
srun: error: hpc-sm-01-15: task 0: Exited with exit code 1
```

### ZC
2026-11-24 14:01:33 rc=0 (srun) ended: Terminating on cycle limit wall=104 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 6.305118e+01 | faile
- deck vs production: > final_restart = false;> nlim = 200;> dt_report = 0;> vet_gd_gl_nphi = 36;> dt = 1e30;> dt = 1e30;> dt = 1e30;> dt = 1e30;> dt = 1e30;
- first lines: <rad_vet> vet_gd direction set: 252 directions; moment errors |sum w - 1| 3.33067e-15, |sum w n| 3.52908e-17, |sum w nn - I/3| 2.77556e-16 (rotated 1.11022e-16), 4th-moment anisotropy 5.82867e-16 (rot
```
/resnick/groups/carnegie_poc/jingze/sdo_1124/ZC/bin/m35zams.vet.00000.bin  t=0.000000e+00  shells 384, thin 76 (i 308-383)
  F3/|F1| odd: thin max 1.050e-02 median 6.615e-03 | thick(i<thin) max 4.797e-04 | thin mean <p> 2.43e-10
  v_phi rms: thin 0.0000e+00  mass-weighted all 0.0000e+00  max|v_phi| thin 0.0000e+00
/resnick/groups/carnegie_poc/jingze/sdo_1124/ZC/bin/m35zams.vet.00001.bin  t=6.961010e+02  shells 384, thin 76 (i 308-383)
  F3/|F1| odd: thin max 9.169e-03 median 5.897e-03 | thick(i<thin) max 1.260e-04 | thin mean <p> 4.01e-05
  v_phi rms: thin 1.5540e+04  mass-weighted all 4.7043e+00  max|v_phi| thin 1.4053e+05
Traceback (most recent call last):
  File "/resnick/groups/carnegie_poc/jingze/sdo_1124/raymeas.py", line 28, in <module>
    f1, f3, e = full(d, 'vet_f1'), full(d, 'vet_f3'), full(d, 'vet_e')
  File "/resnick/groups/carnegie_poc/jingze/sdo_1124/raymeas.py", line 15, in full
    x2 = np.array(d['mb_x2f']); x3 = np.array(d['mb_x3f']); a = d['mb_data'][n]
KeyError: 'vet_f1'
srun: error: hpc-sm-01-15: task 0: Exited with exit code 1
Traceback (most recent call last):
  File "/resnick/groups/carnegie_poc/jingze/sdo_1124/raymeas.py", line 28, in <module>
    f1, f3, e = full(d, 'vet_f1'), full(d, 'vet_f3'), full(d, 'vet_e')
  File "/resnick/groups/carnegie_poc/jingze/sdo_1124/raymeas.py", line 15, in full
    x2 = np.array(d['mb_x2f']); x3 = np.array(d['mb_x3f']); a = d['mb_data'][n]
KeyError: 'vet_f1'
srun: error: hpc-sm-01-15: task 0: Exited with exit code 1
```

## S1 (bsg, bsg3d.pre-compton.rst) (11-24 14:03)
| arm | rows | s/cyc (cycles) | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | L_top/L_in end | L_bot/L_in end | KE | dE/E | dM/M | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries | SDO: lines, K1K2<0 mean share, max ratio, cells>0.1, >1, ef<=0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S1off | 300 | 0.6436 (c6950-7200) | 146.7 / 152.8 / 164.7 | 15.0, 3.39 | 14.3, 3.92 | 1.15812 | 1.01394 | 1.31656e+42 | -1.53e-05 | -2.95e-15 | 0 / 0 / 0 / 0 | 0 / 0 / 0 | 300, 0.5943, 9.863e-05, 0, 0, 0 |
| S1on | 300 | 0.6427 (c6950-7200) | 146.7 / 152.8 / 164.7 | 14.9, 3.39 | 14.3, 3.92 | 1.15813 | 1.01395 | 1.31667e+42 | -1.53e-05 | -3.10e-16 | 0 / 0 / 0 / 0 | 0 / 0 / 0 | 300, 0.5941, 4.047e-01, 11, 0, 0 |

- last common row 300: L_top 4.78e-06 L_in 0.00e+00 KE 8.45e-05 Etot -2.84e-10 M 2.33e-15 Lz -1.71e-04 (on - off relative)
- **PASS** (rule: extra CW/fallback/BE/retry, ef<=0, dt mean < 0.95 off, s per sim s > 1.05 off, |L_top| > 1 %) 
- ends: S1off: 2026-11-24 13:59:52 rc=0 wall=222 s CW 0 FATAL 0 inf 0 DT_RETRY 0 ; S1on: 2026-11-24 14:03:40 rc=0 wall=228 s CW 0 FATAL 0 inf 0 DT_RETRY 0

## zzzzh items (11-24 14:04)
- **S4 dropped:** both S4 arms stopped at ~13:57-14:04 (STOP file, tests); not scored. Partial dirs are kept without rst.
- **ZA's 'cycle -21965':** an artefact of my live table (it subtracts the S2 base cycle 22160 from every t_Z row). ZA is a fresh start: first hst row t = 0.0, dt 3.4805, cycle 0; it ran its first 5 cycles and ended on the cycle limit (200).
- Z arms ran ZA, ZB, ZC one after the other on 01-15 (already running when the node request came): all ended on the cycle limit, srun exit 0, CW 0, FATAL 0, DT_RETRY 0; ZB (S10, vet_gd_ls = 10, GL keys removed) did NOT FATAL.
- Ray effect |F3/F1| odd, thin shells, vet dump (raymeas.py, i 308-383): ZA 7x18 t=0 max 1.485e-2 median 1.13e-2, t=696 s 1.385e-2 / 1.065e-2; ZB S10 2.504e-2 / 2.48e-2, then 2.362e-2 / 2.33e-2; ZC 7x36 1.050e-2 / 6.6e-3, then 9.17e-3 / 5.9e-3. v_phi rms in the thin layer at c200 ZA 1.37e4, ZB 2.86e4, ZC 1.55e4 cm/s. s/cyc 0.340 / 0.341 / 0.343; dt identical; L_top/L_in 1.00046 in all. (The KeyError lines in rays.txt are my script also feeding the vet_vet dump; no data loss.)
- **S1 (bsg) PASS** (0.6436 / 0.6427 s/cyc, dt 152.8 identical, L_top +4.8e-6, SDO max ratio off 1e-4, on 0.40, 11 cells > 0.1, none > 1, ef<=0 0). S2, S3, S5, S1 all PASS; S2L is running (off on 01-15, on on 02-16 after S4 ended).

## S2L (M35ZAMS, 900 cycles; off on 01-15, on on 02-16) (11-24 14:10)
| arm | rows | s/cyc (cycles) | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | L_top/L_in end | L_bot/L_in end | KE | dE/E | dM/M | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries | SDO: lines, K1K2<0 mean share, max ratio, cells>0.1, >1, ef<=0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S2Loff | 900 | 0.2941 (c22160-23050) | 3.481 / 3.481 / 3.481 | 18.4, 2.00 | 15.5, 2.21 | 1.00054 | 1.00001 | 5.56605e+35 | -8.26e-08 | -2.11e-15 | 0 / 0 / 0 / 0 | 0 / 0 / 0 | 900, 0.3173, 1.019e-05, 0, 0, 0 |
| S2Lon | 900 | 0.2925 (c22160-23050) | 3.481 / 3.481 / 3.481 | 18.3, 2.00 | 15.6, 2.21 | 1.00054 | 1.00001 | 5.56705e+35 | -8.26e-08 | 1.74e-15 | 0 / 0 / 0 / 0 | 0 / 0 / 0 | 900, 0.3171, 4.676e-01, 497, 0, 0 |

- last common row 900: L_top -6.83e-07 L_in 0.00e+00 KE 1.79e-04 Etot -2.97e-11 M 9.92e-16 Lz 1.03e-03 (on - off relative)
- **PASS** (rule: extra CW/fallback/BE/retry, ef<=0, dt mean < 0.95 off, s per sim s > 1.05 off, |L_top| > 1 %) 
- ends: S2Loff: 2026-11-24 14:08:11 end: t_S2Loff|COMPLETED|0:0|00:04:44; log: Terminating on cycle limit CW 0 FATAL 0 DT_RETRY 0 ; S2Lon: 2026-11-24 14:10:52 rc=0 wall=287 s CW 0 FATAL 0 inf 0 DT_RETRY 0

### S2 / S2L log ends (zzzzf): last 3 lines of log.0 and Slurm exit status
- **S2off**: t_S2off|COMPLETED|0:0|00:01:53
```
<rad_vet> time_scheme=hesdirk2: stage steps=3.000000e+02 backward-Euler steps=0.000000e+00 stage fallbacks=0.000000e+00 vet clips=0.000000e+00 face_fct calls=6.000000e+02 skipped(no face flux)=0.00000
<rad_vet> time2_auto_be=16: skipped-stage BE steps=0.000000e+00 stage retries=0.000000e+00 (succeeded 0.000000e+00)
<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 8.812223e+01 | failed stage attempts 0.000000e+00 | their BE redo 0.000000e+00 | auto_be steps 0.000000e+00 | other BE steps 0.000000e+00 | 
```
- **S2on**: t_S2on|COMPLETED|0:0|00:01:52
```
<rad_vet> time_scheme=hesdirk2: stage steps=3.000000e+02 backward-Euler steps=0.000000e+00 stage fallbacks=0.000000e+00 vet clips=0.000000e+00 face_fct calls=6.000000e+02 skipped(no face flux)=0.00000
<rad_vet> time2_auto_be=16: skipped-stage BE steps=0.000000e+00 stage retries=0.000000e+00 (succeeded 0.000000e+00)
<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 8.863439e+01 | failed stage attempts 0.000000e+00 | their BE redo 0.000000e+00 | auto_be steps 0.000000e+00 | other BE steps 0.000000e+00 | 
```
- **S2Loff**: t_S2Loff|COMPLETED|0:0|00:04:44
```
<rad_vet> time_scheme=hesdirk2: stage steps=9.000000e+02 backward-Euler steps=0.000000e+00 stage fallbacks=0.000000e+00 vet clips=0.000000e+00 face_fct calls=1.800000e+03 skipped(no face flux)=0.00000
<rad_vet> time2_auto_be=16: skipped-stage BE steps=0.000000e+00 stage retries=0.000000e+00 (succeeded 0.000000e+00)
<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 2.543673e+02 | failed stage attempts 0.000000e+00 | their BE redo 0.000000e+00 | auto_be steps 0.000000e+00 | other BE steps 0.000000e+00 | 
```
- **S2Lon**: t_S2Lon|COMPLETED|0:0|00:04:46
```
<rad_vet> time_scheme=hesdirk2: stage steps=9.000000e+02 backward-Euler steps=0.000000e+00 stage fallbacks=0.000000e+00 vet clips=0.000000e+00 face_fct calls=1.800000e+03 skipped(no face flux)=0.00000
<rad_vet> time2_auto_be=16: skipped-stage BE steps=0.000000e+00 stage retries=0.000000e+00 (succeeded 0.000000e+00)
<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 2.530778e+02 | failed stage attempts 0.000000e+00 | their BE redo 0.000000e+00 | auto_be steps 0.000000e+00 | other BE steps 0.000000e+00 | 
```
max|ef/es-1| per 50 cycles (S2Lon): 0.12 0.14 0.15 0.24 0.29 0.47 0.38 0.34 0.36 0.30 0.32 0.25 0.19 0.32 0.40 0.28 0.33 0.47 

## Marked per your 24zzzzj: all S passes are CONDITIONAL
- **S1, S2, S2L, S3, S5: 'PASS' = conditional: old binary (a840d7e92), ef<=0 NOT tested.** The key has a positivity defect in floorcap gas (the increment is built before FloorGasCap) that the ef<=0 counter does not see, so 'ef<=0 0' in these tables is not a valid test. The other criteria (CW, fallback, BE, retry, dt, cost, L_top) stand as measured. S2L ran to its end (posted above); no further S arm; S3 and S5 are rerun only when you give the fixed binary.

## S3F (fixed binary sdorder2 09d1e4df4, both arms, SDO_PRINT=1), 11-24 14:58
| arm | rows | s/cyc (cycles) | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | L_top/L_in end | L_bot/L_in end | KE | dE/E | dM/M | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries | SDO: lines, K1K2<0 mean share, max ratio, cells>0.1, >1, ef<=0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S3Foff | 300 | 1.0333 (c7000-7250) | 9.862 / 11.5 / 12.01 | 34.1, 4.95 | 16.1, 4.01 | 0.99958 | 1.00016 | 1.57167e+39 | 9.32e-07 | -3.78e-16 | 0 / 0 / 0 / 0 | 0 / 0 / 0 | 300, 0.3227, 1.882e-02, 0, 0, 0 |
| S3Fon | 300 | 1.0325 (c7000-7250) | 9.86 / 11.5 / 12.01 | 34.3, 4.92 | 16.0, 3.97 | 0.99975 | 1.00016 | 1.57194e+39 | 9.06e-07 | -1.89e-15 | 0 / 0 / 0 / 0 | 0 / 0 / 0 | 300, 0.3234, 2.277e+00, 95623, 255, 0 |

- last common row 300: L_top 1.74e-04 L_in 6.91e-16 KE 1.71e-04 Etot -2.59e-08 M 1.89e-15 Lz 1.18e-05 (on - off relative)
- floor_pay sum (CLEAN FALLBACK lines): off 1.781e+26 on 1.719e+26 ratio 0.965; per 50 cycles off: 4.51e+25 3.56e+25 4.23e+25 3.15e+25 1.59e+25 6.86e+24 7.92e+23 | on: 4.57e+25 3.44e+25 4.17e+25 3.23e+25 1.32e+25 4.45e+24 1.44e+23
- **PASS** (rule: extra CW/fallback/BE/retry, ef<=0, dt mean < 0.95 off, s per sim s > 1.05 off, |L_top| > 1 %) 
- ends: S3Foff: 2026-11-24 14:52:20 rc=0 wall=346 s CW 0 FATAL 0 inf 0 DT_RETRY 0 | S3Fon: 2026-11-24 14:58:01 rc=0 wall=341 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- ef<=0 per 50 cycles (post-cap) off/on: 0 0 0 0 0 0 ; 0 0 0 0 0 0 ; 

## S5F (fixed binary sdorder2 09d1e4df4, both arms, SDO_PRINT=1), 11-24 15:16
| arm | rows | s/cyc (cycles) | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | L_top/L_in end | L_bot/L_in end | KE | dE/E | dM/M | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries | SDO: lines, K1K2<0 mean share, max ratio, cells>0.1, >1, ef<=0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S5Foff | 300 | 3.0398 (c5300-5550) | 10.51 / 32.88 / 46.91 | 14.2, 4.84 | 8.4, 3.61 | 2.44083 | 0.98459 | 1.78357e+44 | -6.01e-05 | 3.09e-16 | 0 / 0 / 0 / 0 | 0 / 0 / 0 | 300, 0.2879, 9.512e+02, 9088084, 502823, 502823 |
| S5Fon | 300 | 3.0016 (c5300-5550) | 12.04 / 38.71 / 50.52 | 14.9, 4.71 | 8.7, 3.34 | 2.42044 | 0.98680 | 1.78883e+44 | -7.02e-05 | 1.55e-16 | 0 / 0 / 0 / 0 | 0 / 0 / 0 | 300, 0.2876, 7.910e+01, 27727803, 13085411, 0 |

- last common row 300: L_top -8.35e-03 L_in 1.72e-16 KE 2.95e-03 Etot -1.02e-05 M -1.55e-15 Lz -7.91e-07 (on - off relative)
- floor_pay sum (CLEAN FALLBACK lines): off 3.280e+32 on 5.000e+27 ratio 0.000; per 50 cycles off: 2.04e+32 1.24e+32 8.71e+26 8.57e+26 8.66e+26 8.01e+26 4.23e+26 | on: 1.44e+27 2.31e+26 5.96e+26 6.88e+26 8.15e+26 7.64e+26 4.63e+26
- **PASS** (rule: extra CW/fallback/BE/retry, ef<=0, dt mean < 0.95 off, s per sim s > 1.05 off, |L_top| > 1 %) 
- ends: S5Foff: 2026-11-24 15:16:21 rc=0 wall=977 s CW 0 FATAL 0 inf 0 DT_RETRY 0 | S5Fon: 2026-11-24 15:16:08 rc=0 wall=964 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- ef<=0 per 50 cycles (post-cap) off/on: 502890 1 1 1 0 0 ; 0 0 0 0 0 0 ; 

## S5F matched-time remark
- S5F on vs off at matched time (300 rows): L_top rel end 4.504e-05 max 2.270e-04
- S5F on vs off at matched time (300 rows): L_bot rel end -5.123e-07 max 1.182e-05
- S5F on vs off at matched time (300 rows): KE_int rel end -3.420e-05 max 3.420e-05
- model time advanced: off 9848.7 s on 11598.1 s
- The 'last common row' line of the S5F table compares the 300th rows, which sit at different model times (dt differs by 18 %), so its -8.4e-3 in L_top mixes evolution; the lines above are the matched-time numbers. The off arm's ef<=0 502823 / max ratio 951 are printed with the key off (the print uses the key's state).  Both are stated as measured; the meaning of the post-cap counter in the off arm is yours to judge.

## S3G (fixed binary sdorder3 95813cf85, both arms, SDO_PRINT=1), 11-24 16:12
| arm | rows | s/cyc (cycles) | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | L_top/L_in end | L_bot/L_in end | KE | dE/E | dM/M | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries | SDO: lines, K1K2<0 mean share, max ratio, cells>0.1, >1, ef<=0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S3Goff | 300 | 1.0152 (c7000-7250) | 9.862 / 11.5 / 12.01 | 34.1, 4.95 | 16.1, 4.01 | 0.99958 | 1.00016 | 1.57167e+39 | 9.32e-07 | -3.78e-16 | 0 / 0 / 0 / 0 | 0 / 0 / 0 | 0, 0.0000, nan, nan, nan, nan |
| S3Gon | 300 | 1.0138 (c7000-7250) | 9.86 / 11.5 / 12.01 | 34.2, 4.96 | 16.1, 3.97 | 0.99977 | 1.00016 | 1.57194e+39 | 9.06e-07 | -2.27e-15 | 0 / 0 / 0 / 0 | 0 / 0 / 0 | 0, 0.0000, nan, nan, nan, nan |

- last common row 300: L_top 1.93e-04 L_in 6.91e-16 KE 1.71e-04 Etot -2.59e-08 M 1.51e-15 Lz -4.25e-06 (on - off relative)
- floor_pay sum (CLEAN FALLBACK lines): off 1.781e+26 on 1.708e+26 ratio 0.959; per 50 cycles off: 4.51e+25 3.56e+25 4.23e+25 3.15e+25 1.59e+25 6.86e+24 7.92e+23 | on: 4.55e+25 3.43e+25 3.95e+25 3.30e+25 1.34e+25 4.79e+24 2.22e+23
- **PASS** (rule: extra CW/fallback/BE/retry, ef<=0, dt mean < 0.95 off, s per sim s > 1.05 off, |L_top| > 1 %) 
- ends: S3Goff: 2026-11-24 16:06:53 rc=0 wall=336 s CW 0 FATAL 0 inf 0 DT_RETRY 0 | S3Gon: 2026-11-24 16:12:29 rc=0 wall=336 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- ef<=0 per 50 cycles (post-cap) off/on: 0 ; 0 ; 

## S5G (fixed binary sdorder3 95813cf85, both arms, SDO_PRINT=1), 11-24 16:17
| arm | rows | s/cyc (cycles) | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | L_top/L_in end | L_bot/L_in end | KE | dE/E | dM/M | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries | SDO: lines, K1K2<0 mean share, max ratio, cells>0.1, >1, ef<=0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S5Goff | 300 | 2.9995 (c5300-5550) | 10.51 / 32.88 / 46.91 | 14.2, 4.84 | 8.4, 3.61 | 2.44083 | 0.98459 | 1.78357e+44 | -6.01e-05 | 3.09e-16 | 0 / 0 / 0 / 0 | 0 / 0 / 0 | 0, 0.0000, nan, nan, nan, nan |
| S5Gon | 300 | 2.9599 (c5300-5550) | 11.63 / 37.35 / 49.75 | 14.7, 4.74 | 8.5, 3.39 | 2.42542 | 0.98681 | 1.78769e+44 | -6.79e-05 | 2.63e-15 | 0 / 0 / 0 / 0 | 0 / 0 / 0 | 0, 0.0000, nan, nan, nan, nan |

- last common row 300: L_top -6.31e-03 L_in 0.00e+00 KE 2.31e-03 Etot -7.82e-06 M 9.27e-16 Lz -5.80e-07 (on - off relative)
- floor_pay sum (CLEAN FALLBACK lines): off 3.280e+32 on 9.917e+32 ratio 3.023; per 50 cycles off: 2.04e+32 1.24e+32 8.71e+26 8.57e+26 8.66e+26 8.01e+26 4.23e+26 | on: 5.78e+26 9.92e+32 2.07e+25 2.06e+25 1.98e+25 1.72e+25 9.90e+24
- **FAIL: floor_pay on 9.917e+32 > 1.2 x off 3.280e+32** (rule: extra CW/fallback/BE/retry, ef<=0, dt mean < 0.95 off, s per sim s > 1.05 off, |L_top| > 1 %) 
- ends: S5Goff: 2026-11-24 16:17:21 rc=0 wall=964 s CW 0 FATAL 0 inf 0 DT_RETRY 0 | S5Gon: 2026-11-24 16:17:13 rc=0 wall=956 s CW 0 FATAL 0 inf 0 DT_RETRY 0
- ef<=0 per 50 cycles (post-cap) off/on: 0 ; 0 ; 
