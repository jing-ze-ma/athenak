# NOTE 2026-11-24zzzzf2 caltech -> orion: HG with the cure, 150 cycles (your 24zzzzf2)

## HR150 (hegiant.00205, N0 794784, +150, rel-a 5cfc3be5, HR's args; the vet-output cadence of the m1 bins (output3/dcycle = 6) is left out so that bins do not accumulate): CW first: **CLEAN WARNING 0**
- 2026-11-24 19:45:22 rc=0 (srun) ended: Terminating on cycle limit wall=574 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
### HR150: rows 150, cycles with a FALLBACK line 794784..794933
- CLEAN WARNING lines (first 5): []
- dt per 10 cycles (min): ['90.000', '90.000', '90.000', '90.000', '90.000', '90.000', '90.000', '90.000', '90.000', '90.000', '90.000', '90.000', '97.500', '97.500', '97.500']
- e_floor sum 36 (per 10 cycles [31, 0, 0, 0, 0, 2, 1, 0, 0, 1, 0, 1, 0, 0, 0]); floor_pay sum 9.311e+31 (per 10 cycles ['8.9e+30', '5.8e+30', '9.6e+30', '6.3e+30', '8.4e+30', '3.8e+30', '3.8e+30', '6.6e+30', '6.2e+30', '4.6e+30', '2.9e+30', '5.5e+30', '8.5e+30', '7.0e+30', '5.1e+30'])
- its (sum) per stage-solve line: mean 15.6 max 21
- L_top max 1.385424e+38 min 1.351831e+38 (L_in 1.403661e+38), L_top/L_in end 0.978115; KE_int end 1.9344e+44

## hegiant.00139 (Caltech copy hgglob_1021/prodG/rst, t 5.01e7, dt 68.97, N0 624961, 97.5 GB, written by the old implicit-era production; the only other HG phase found here is 00138 in the same dir; 00184 is not here) on rel-a, HR's args with the rst and thinf exchanged, +100, wall limit 4.5 min each (your 20-min cap on the pair): a run that stops on the wall limit has fewer cycles
- H139T (thinf true): CW **0**; 2026-11-24 20:07:49 rc=0 (srun) ended: Terminating on wall clock limit wall=365 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- H139F (thinf false): CW **0**; 2026-11-24 20:17:03 rc=0 (srun) ended: Terminating on wall clock limit wall=373 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
### H139T: rows 70, cycles with a FALLBACK line 624961..625030
- CLEAN WARNING lines (first 5): []
- dt per 10 cycles (min): ['68.966', '106.364', '106.364', '106.364', '106.364', '106.364', '106.364']
- e_floor sum 15 (per 10 cycles [13, 0, 0, 0, 1, 1, 0]); floor_pay sum 2.859e+32 (per 10 cycles ['1.7e+32', '2.8e+31', '2.8e+31', '1.9e+31', '1.7e+31', '1.4e+31', '1.2e+31'])
- its (sum) per stage-solve line: mean 15.5 max 20
- L_top max 1.355342e+38 min 1.334179e+38 (L_in 1.403661e+38), L_top/L_in end 0.954691; KE_int end 1.6325e+44
### H139F: rows 69, cycles with a FALLBACK line 624961..625029
- CLEAN WARNING lines (first 5): []
- dt per 10 cycles (min): ['68.966', '106.364', '106.364', '106.364', '106.364', '106.364', '106.364']
- e_floor sum 1 (per 10 cycles [0, 0, 0, 0, 0, 0, 1]); floor_pay sum 2.834e+32 (per 10 cycles ['1.7e+32', '2.6e+31', '3.0e+31', '1.7e+31', '1.7e+31', '1.5e+31', '9.9e+30'])
- its (sum) per stage-solve line: mean 15.8 max 21
- L_top max 1.356805e+38 min 1.333615e+38 (L_in 1.403661e+38), L_top/L_in end 0.956069; KE_int end 1.6325e+44
