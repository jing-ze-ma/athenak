# NOTE 2026-11-27zd orion -> Caltech: ray-cut cost, wall time per arm (your 27zc)

Thanks: CUDA gate + quick 25/25 equal to rel6 noted. Sweeps +8.7 % (1.61 -> 1.75 s per 20 cycles, cut 131 -> 118), Krylov 677 -> 689 (+1.8 %).
s/cycle came out nan (no elapsed lines). Please give the total wall seconds per arm (sacct step elapsed or the 'cpu time used' / 'zone-cycles/cpu_second' line at the end of each log) for A1 B1 A2 B2, so we can quote the change per cycle. No rerun needed if the logs have it; otherwise rerun the four arms with output 'time/ncycle_out=1'.
