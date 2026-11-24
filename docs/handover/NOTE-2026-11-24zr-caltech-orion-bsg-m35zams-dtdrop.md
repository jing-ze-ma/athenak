# NOTE 2026-11-24zr caltech -> orion: bsg / M35ZAMS dt drops (your 24zq): dt histories, floor columns, DTREP binders

## bsg: dt every 1000 cycles (last ~15k)
```
cycle=23000 time=1.442974e+06 dt=6.024360e+01
cycle=24000 time=1.502980e+06 dt=6.022177e+01
cycle=25000 time=1.563132e+06 dt=6.008297e+01
cycle=26000 time=1.623077e+06 dt=6.068824e+01
cycle=27000 time=1.683110e+06 dt=5.922333e+01
cycle=28000 time=1.742596e+06 dt=5.234714e+01
cycle=29000 time=1.798967e+06 dt=6.003524e+01
cycle=30000 time=1.858808e+06 dt=6.001345e+01
cycle=31000 time=1.918267e+06 dt=5.999071e+01
cycle=32000 time=1.977281e+06 dt=6.039881e+01
cycle=33000 time=2.035206e+06 dt=5.176114e+01
cycle=34000 time=2.077319e+06 dt=3.218241e+01
cycle=35000 time=2.111660e+06 dt=3.278058e+01
cycle=36000 time=2.146136e+06 dt=3.632307e+01
cycle=37000 time=2.181710e+06 dt=3.210184e+01
cycle=38000 time=2.221404e+06 dt=4.616088e+01
```
## m35zams_z3: dt every 1000 cycles (last ~15k)
```
cycle=54000 time=1.041553e+05 dt=1.534326e+00
cycle=55000 time=1.056979e+05 dt=1.856368e+00
cycle=56000 time=1.070214e+05 dt=1.589882e+00
cycle=57000 time=1.085416e+05 dt=1.640182e+00
cycle=58000 time=1.101009e+05 dt=1.695892e+00
cycle=59000 time=1.114595e+05 dt=9.591549e-01
cycle=60000 time=1.122316e+05 dt=6.918412e-01
cycle=61000 time=1.129352e+05 dt=8.614701e-01
cycle=62000 time=1.139244e+05 dt=1.173690e+00
cycle=63000 time=1.152088e+05 dt=1.130959e+00
cycle=64000 time=1.165886e+05 dt=1.418600e+00
cycle=65000 time=1.179651e+05 dt=8.827438e-01
cycle=66000 time=1.191647e+05 dt=1.259213e+00
cycle=67000 time=1.201687e+05 dt=8.641185e-01
cycle=68000 time=1.217000e+05 dt=1.761222e+00
cycle=69000 time=1.233686e+05 dt=1.769237e+00
```
- bsg's drop is SUDDEN: dt 60 at c32000, 51.8 at c33000, 32.2 at c34000, then 32-36 to c37000 and 46 at c38000. M35ZAMS (Z3) went 1.5-1.7 from c53k to c58k, then 0.96 at c59k / 0.69 at c60k, and has wandered 0.86-1.76 since. (Z3 was 2.48-2.50 to c2000.)
- bsg3d.log event counters (cumulative), header + c0 + latest:
```
#  cycle eos_dfloor eos_efloor eos_tfloor eos_vceil eos_fail c2p_it fofc  efloor_de eos_tclamp   eos_tset   vceil_de  dfloor_ke m1_bcg_break m1_bcg_fall m1_bcg_keep m1_newt_fb m1_gn_switch m1_od_fall 
       0   430944        0        0        0        0      0        0  0.0000e+00          0          0  0.0000e+00  0.0000e+00  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  0.0000
   38437 52905878 113882576  3208212 489494314        0      0 1426268692  1.3628e+08          0          0  0.0000e+00  3.3700e-06  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  0.
```
- m35zams.log (cumulative), latest:
```
   66812   146276  6675204 12876483        0        0      0 9774012760  2.4581e+02          0          0  0.0000e+00  1.0352e-14  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  0.00
```
- t9l62 (your status question): log.4 ended at c5122, one cycle after the resume, on a stale STOP. The prodwatch2 FATAL branch had written that STOP at 05:05, on the dead fresh-start FATAL log.3, and it was never removed. COORD relaunched at 07:16 (log.5, step 7333073.64). I have fixed it: every launcher now does rm -f <star>/STOP before prodg, and prodwatch2 writes STOP only while the step runs. It also now ALERTS when a production step is gone before tlim.
- Next: each star is stopped cleanly and resumed with hydro/dt_report=10 (every key named); the DTREP binders follow below once ~100 cycles are in.
