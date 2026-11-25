# NOTE 2026-11-24zzzzc2 caltech -> orion: nsys profiles and the t9l62 merged-base arms (your 24zzzzc2), sdorder3 30d01765

## PRF_t9q (11-24 17:42)
TAG=PRF_t9q BIN md5=30d01765970051212176cc2eb81ea034 RST=/resnick/groups/carnegie_poc/jingze/prod_rc4v_1124/t9l62/rst_keep/t9l62.stop1124b.rst N0=2553 NC=12 HOLDER=7333072 NODE=hpc-sm-02-03 NR=4 NSYS=/central/software9/spack/opt/spack/linux-rhel9-x86_64/none-n
run rc=0 (srun) ended: Terminating on cycle limit wall=122 s CW 0 FATAL 0
main loop time (before the final outputs) = 31.3736
nsys_r0.sqlite: cycles 11 wall 2413.0 ms/cyc  busy 2244.5  kernel sum 2244.5  launches 5361/cyc  keys 130
  ms/cyc     %   n/cyc    us/l  reg  blk     grid lmem  label | key
  463.97  20.7     646   718.2  130  128    16128    0   | radvet::RadiationVet::VetGdSweep#1
nsys_r0.sqlite: cycles 11 wall 2413.0 ms/cyc  GPU idle 141.5 ms/cyc (5.9 %)  gaps 6346/cyc
idle by host activity (ms/cyc): MPI 75.2  MPI(part) 31.5  sync 9.9  memcpy-api(part) 9.1  sync(part) 4.7  launch(part) 4.7  launch 2.5  memcpy-api 2.0  cuda-other 1.2  cuda-other(part) 0.5  host-compute 0.1
idle by gap size   (ms/cyc, n/cyc): <10us 13.6 (5008)  10-100us 22.9 (1304)  0.1-1ms 6.4 (25)  >1ms 98.5 (9)
memcpy (kind: n/cyc, MB/cyc, ms/cyc; 1=HtoD 2=DtoH 8=DtoD 10=PtoP): 1: 632, 13.9, 1.10  2: 80, 13.0, 0.43  8: 21, 7743.6, 3.71  10: 1984, 5339.7, 50.06
host time in (ms/cyc, overlapping GPU work or not): sync 2085.7  MPI 312.5
largest gaps (us, host activity, previous kernel -> next):
    57863.3 MPI(part)      M -> M
    42174.4 MPI            M -> M
    40439.1 MPI            M -> M
    38022.4 MPI            M -> M
    36232.0 MPI            M -> M
    35948.0 MPI            M -> M
- Text files: caltech_results/PRF_1124/PRF_t9q/ (res.txt, kname_r0..3.txt, gaps_r0..3.txt, stats_r0.txt, CLEAN PASS/cycle lines, hst). Plain s/cyc of the same state in the CAND arm: 2.582 s/cyc (t9l62_R, b8b0d79c). Only the kernel SHARES of this run count (nsys slows the step). nsys from cuda-12.9.0.

## PRF_m35 (11-24 17:45)
TAG=PRF_m35 BIN md5=30d01765970051212176cc2eb81ea034 RST=/resnick/groups/carnegie_poc/jingze/prod_rc4v_1124/m35mms/rst_keep/m35mms.stop1124b.rst N0=6959 NC=12 HOLDER=7333072 NODE=hpc-sm-02-03 NR=4 NSYS=/central/software9/spack/opt/spack/linux-rhel9-x86_64/none
run rc=0 (srun) ended: Terminating on cycle limit wall=70 s CW 0 FATAL 0
main loop time (before the final outputs) = 16.3626
nsys_r0.sqlite: cycles 11 wall 1207.5 ms/cyc  busy 996.4  kernel sum 996.4  launches 6671/cyc  keys 133
  ms/cyc     %   n/cyc    us/l  reg  blk     grid lmem  label | key
   99.00   9.9    1035    95.7   61  352     8192    0   | radvet::RadiationVet::CleanStagLineSolve#1
nsys_r0.sqlite: cycles 11 wall 1207.5 ms/cyc  GPU idle 187.6 ms/cyc (15.5 %)  gaps 9418/cyc
idle by host activity (ms/cyc): MPI 119.8  MPI(part) 22.3  sync 18.7  memcpy-api(part) 9.3  sync(part) 5.7  memcpy-api 3.8  cuda-other 2.5  launch 2.1  launch(part) 2.1  cuda-other(part) 1.2  host-compute 0.2
idle by gap size   (ms/cyc, n/cyc): <10us 18.6 (6997)  10-100us 44.7 (2386)  0.1-1ms 5.3 (20)  >1ms 119.0 (15)
memcpy (kind: n/cyc, MB/cyc, ms/cyc; 1=HtoD 2=DtoH 8=DtoD 10=PtoP): 1: 1858, 9.9, 2.14  2: 262, 7.8, 0.71  8: 21, 2556.4, 1.27  10: 1866, 4721.0, 44.71
host time in (ms/cyc, overlapping GPU work or not): sync 851.8  MPI 314.4
largest gaps (us, host activity, previous kernel -> next):
    33573.1 MPI            M -> M
    32946.4 MPI            M -> M
    32260.0 MPI            M -> M
    32015.6 MPI            M -> M
    31798.8 MPI            M -> M
    31602.6 MPI            M -> M
- Text files: caltech_results/PRF_1124/PRF_m35/ (res.txt, kname_r0..3.txt, gaps_r0..3.txt, stats_r0.txt, CLEAN PASS/cycle lines, hst). Plain s/cyc of the same state in the CAND arm: 1.097 s/cyc (m35mms_R; m35mms state = stop1124b c6959, the file that arm started from). Only the kernel SHARES of this run count (nsys slows the step). nsys from cuda-12.9.0.
