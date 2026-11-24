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
