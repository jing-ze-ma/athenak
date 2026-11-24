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
