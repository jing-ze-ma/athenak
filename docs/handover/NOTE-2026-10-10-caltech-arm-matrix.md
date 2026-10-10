# NOTE Caltech -> viper: xthinfix-1009 ARM MATRIX (413b38af), 10-10

Answers TASK-2026-10-10-caltech-arm-matrix.md. Run trees stay on Caltech: /resnick/groups/carnegie_poc/jingze/arms_1010
(run/, logs/; REPO = git archive of 413b38af in repo/).

STATUS: complete. Jobs:
- GPU: 4309030 (smoke 4308662).
- CPU: order 10 x 32 cores 4308649-58 (1710 runs, all rc 0), gates 4308659, beams 4308660 (50 cores), ana 4308661.
- Debug build + gdb: 4309595.

## Verdict

1. Order (1-D): sq = sqf = sqv = sqa (bitwise in 1-D).
   - sq matches central on the thin pulse (k0.128 hesdirk2 X 1.94 2.04) and rw_t10 X converges (1.97 1.88).
   - Hopf L1 at 512 is 3.87e-4 (hrx0 3.71e-4).
   - atm_self converges in every region: tau<0.1 1.12, tau 1-10 1.15, thick 1.60, Kn>0.3 0.95 (confirms viper).
   - aq / knp match sq on the atmosphere (3.88e-4, Kn>0.3 1.00), but the thin pulse breaks in X (0.69; C 1.72); aq also fails rw_t10 (X -1.00, C -0.47).
   - kn01 / kn03 / kn1 (dc) do not converge in Kn>0.3 (orders ~0 / ~0 / 0.47).
2. Gates: EVERY plm arm (sq, aq, knp) aborts in G1 / G3 (2-D, nx2 4): 12 of 12 runs give rc 139. hr and kn03 pass G1 / G3 / G5, and G5 passes for all arms.
   - Root cause (Debug build + gdb, same on G1 and G3): ImplicitMusclBuild (rad_m1_implicit.cpp:12490) calls ImplicitHaloExchange(3, ...). nq = 3 is not one of the widths there (M1_NHALO_T 14, M1_NHALO_Q 4, M1_NVIMP_X 12), so it falls to the else branch: the width-1 krw scratch. ImplicitHaloCopy (:2741) then writes component n = 1 into a view of extent [1,1,1,8,132] ("Kokkos::View ERROR: out of bounds access ... indices [0,1,0,0,0] but extents [1,1,1,8,132]").
   - In Release this is a segfault. It hits only the generic CC-boundary halo path; the GPU runs use the implicit_halo_mpi path and do not crash.
3. Beams (no FATAL in hr / aq / kn* / knp; NON-CONVERGED 0 wherever a run finished):
   - sq: xb20 Picard 65 (max 177), cyl 9.8, shd3b 4.7; ba0 DIVERGES (cycle 4) and ba20 DIVERGES (cycle 88). The ba0/ba20 divergence is new versus viper.
   - sqv behaves the same as sq (ba0 diverges at cycle 4, ba20 at cycle 75).
   - sqa diverges on xb20 / ba0 / ba20. sqf diverges on xb20 / ba0 / ba20 / cyl.
   - aq Picard 7.8 / 15.4 / 22.4 / 33.8 (xb20 / ba0 / ba20 / cyl); knp similar (8.0 / 14.3 / 20.6 / 34.0).
   - kn01 / kn03 / kn1: Picard 4.7-7.1, but inner BiCGStab 55-143 (max 200) versus hr 17-45.
4. GPU AG Car A, 2 H200: base (5b304cf9) == all in bin data; base repeat bitwise.
   - kn03 costs the same as all (0.467 vs 0.465 s/cycle, Picard 4.27 vs 4.33) and differs only near R_ph (at most 2e-2 in F2/F3 at 0.95-1.05).
   - sq / aq / knp cost 3.8-4.1x (1.90 / 1.83 / 1.75 s/cycle, Picard 45-47) and give O(1) differences above 1.05 R_ph.
5. Net: sq has the right convergence but diverges on 2 of 5 beams, costs 4x on the GPU, and its plm path segfaults in 2-D on the generic halo path. The dc kn arms are cheap but do not converge in the thin top.

## Caltech adaptations (same bundle issues as 10-09 / steep; copies, not in place)

- Gate inputs (atm2d, pulse2d, marsh2d) still lack `vet_tensor`. `fixkeys.py` adds any missing command-line key to a copy of the input; only vet_tensor was missing.
- /usr/bin/time is replaced by a shell stand-in with the same output.
- beams() writes `<run>/cmd.txt` for cylref/collref.
- python: spack 3.11 + numpy/scipy/matplotlib + py-h5py. ana_atm_self.py needs only numpy + od.py and ran fine (atm_self_regions.txt).
- order ran as 10 round-robin job files (od.py run, 32 cores each). Copy: arms/battery_arms_caltech.sh. GPU: only the #SBATCH header and srun (gpuwrap.sh) changed.
- Debug build for the backtrace: separate tree inc_cpu_built_in_pgens_dbg, CMAKE_BUILD_TYPE=Debug, athena_cpu_built_in_pgens_413b38af_dbg. gdb logs: run/dbg/{g1_sq_1e2,g3_sq_c1}/gdb.log.

## FATAL / NON-CONVERGED

- Beams, Picard DIVERGED (FATAL):
  - ba0: sq c4, sqv c4, sqa c1, sqf c1
  - ba20: sq c88, sqv c75, sqa c6, sqf c2
  - xb20: sqa c6, sqf c1
  - cyl: sqf c85
- Gates, SIGSEGV rc 139: g1_{sq,aq,knp}_{1e2,1e4}, g3_{sq,aq,knp}_{c1,c10}. The 'g?|...' lines at the top of gates.txt are bash echoing the crashed commands.
- No other FATAL. NON-CONVERGED > 0 only in the diverged runs.

## Debug backtrace (g1_sq_1e2; g3_sq_c1 identical)

```

#0  0x00007ffff768ba6c in __pthread_kill_implementation () from /lib64/libc.so.6
#1  0x00007ffff763e686 in raise () from /lib64/libc.so.6
#2  0x00007ffff7628833 in abort () from /lib64/libc.so.6
#3  0x0000000001705f98 in Kokkos::Impl::host_abort (message=0x7fffffff1c50 "Kokkos::View ERROR: out of bounds access label=(\"\") with indices [0,1,0,0,0] but extents [1,1,1,8,132]") at /resnick/home/jingze/ATHENAK/builds/inc_cpu_built_in_pgens_dbg/src/kokkos/
#4  0x0000000000419022 in Kokkos::abort (message=0x7fffffff1c50 "Kokkos::View ERROR: out of bounds access label=(\"\") with indices [0,1,0,0,0] but extents [1,1,1,8,132]") at /resnick/home/jingze/ATHENAK/builds/inc_cpu_built_in_pgens_dbg/src/kokkos/core/src/Ko
#5  0x0000000000465d00 in Kokkos::Impl::view_verify_operator_bounds<Kokkos::HostSpace, Kokkos::View<double*****, Kokkos::LayoutRight, Kokkos::HostSpace>, Kokkos::Impl::ViewMapping<Kokkos::ViewTraits<double*****, Kokkos::LayoutRight, Kokkos::HostSpace>, void>, 
#6  0x00000000011865b1 in Kokkos::View<double*****, Kokkos::LayoutRight, Kokkos::HostSpace>::operator()<int, int, int, int, int> (this=0x7fffffff1ff8) at /resnick/home/jingze/ATHENAK/builds/inc_cpu_built_in_pgens_dbg/src/kokkos/core/src/View/Kokkos_ViewLegacy.
#7  operator() (__closure=0x7fffffff1fd0, m=0, n=1, k=0, j=0) at /resnick/home/jingze/ATHENAK/builds/inc_cpu_built_in_pgens_dbg/src/src/rad_m1/rad_m1_implicit.cpp:2741
#8  0x0000000001287915 in operator() (__closure=0x7fffffff1fb0, idx=@0x7fffffff1f2c: 8) at /resnick/home/jingze/ATHENAK/builds/inc_cpu_built_in_pgens_dbg/src/src/athena.hpp:333
#9  0x00000000013022be in Kokkos::Impl::ParallelFor<par_for<Kokkos::Serial, radm1::RadiationM1::ImplicitHaloCopy(DvceArray5D<double>&, int, int, bool)::<lambda(int, int, int, int)> >(const std::string&, Kokkos::Serial, int const&, int const&, int const&, int c
#10 0x00000000012e6558 in Kokkos::Impl::ParallelFor<par_for<Kokkos::Serial, radm1::RadiationM1::ImplicitHaloCopy(DvceArray5D<double>&, int, int, bool)::<lambda(int, int, int, int)> >(const std::string&, Kokkos::Serial, int const&, int const&, int const&, int c
#11 0x00000000012d16aa in Kokkos::parallel_for<Kokkos::RangePolicy<Kokkos::Serial>, par_for<Kokkos::Serial, radm1::RadiationM1::ImplicitHaloCopy(DvceArray5D<double>&, int, int, bool)::<lambda(int, int, int, int)> >(const std::string&, Kokkos::Serial, int const
#12 0x0000000001287b62 in par_for<Kokkos::Serial, radm1::RadiationM1::ImplicitHaloCopy(DvceArray5D<double>&, int, int, bool)::<lambda(int, int, int, int)> >(const std::string &, Kokkos::Serial, const int &, const int &, const int &, const int &, const int &, c
#13 0x00000000011873db in radm1::RadiationM1::ImplicitHaloCopy (this=0x2437e90, sc=..., nq=3, c0=39, topack=true) at /resnick/home/jingze/ATHENAK/builds/inc_cpu_built_in_pgens_dbg/src/src/rad_m1/rad_m1_implicit.cpp:2730
#14 0x0000000001187fd8 in radm1::RadiationM1::ImplicitHaloExchange (this=0x2437e90, nq=3, c0=39) at /resnick/home/jingze/ATHENAK/builds/inc_cpu_built_in_pgens_dbg/src/src/rad_m1/rad_m1_implicit.cpp:2787
#15 0x0000000001285153 in radm1::RadiationM1::ImplicitMusclBuild (this=0x2437e90) at /resnick/home/jingze/ATHENAK/builds/inc_cpu_built_in_pgens_dbg/src/src/rad_m1/rad_m1_implicit.cpp:12490
#16 0x00000000012637d9 in radm1::RadiationM1::ImplicitSolve (this=0x2437e90, pdrive=0x224b920, stage=1) at /resnick/home/jingze/ATHENAK/builds/inc_cpu_built_in_pgens_dbg/src/src/rad_m1/rad_m1_implicit.cpp:10202
#17 0x0000000001572e99 in TaskList::AddTask<TaskStatus (radm1::RadiationM1::*)(Driver*, int), radm1::RadiationM1>(TaskStatus (radm1::RadiationM1::*)(Driver*, int), radm1::RadiationM1*, TaskID&)::{lambda(Driver*, int)#1}::operator()(Driver*, int) (__closure=0x2
#18 0x000000000157338c in std::__invoke_impl<TaskStatus, TaskList::AddTask<TaskStatus (radm1::RadiationM1::*)(Driver*, int), radm1::RadiationM1>(TaskStatus (radm1::RadiationM1::*)(Driver*, int), radm1::RadiationM1*, TaskID&)::{lambda(Driver*, int)#1}&, Driver*
#19 0x000000000157324d in std::__invoke_r<TaskStatus, TaskList::AddTask<TaskStatus (radm1::RadiationM1::*)(Driver*, int), radm1::RadiationM1>(TaskStatus (radm1::RadiationM1::*)(Driver*, int), radm1::RadiationM1*, TaskID&)::{lambda(Driver*, int)#1}&, Driver*, i
#20 0x0000000001573107 in std::_Function_handler<TaskStatus (Driver*, int), TaskList::AddTask<TaskStatus (radm1::RadiationM1::*)(Driver*, int), radm1::RadiationM1>(TaskStatus (radm1::RadiationM1::*)(Driver*, int), radm1::RadiationM1*, TaskID&)::{lambda(Driver*
#21 0x00000000007475e4 in std::function<TaskStatus (Driver*, int)>::operator()(Driver*, int) const (this=0x2c585f8, __args#0=0x224b920, __args#1=1) at /resnick/software9/spack/opt/spack/linux-rhel9-x86_64/gcc-13.2.0/gcc-13.2.0-w55nxklkmply3rc66ssu3nhqt2d2eg6l/
#22 0x0000000000746fa7 in Task::operator() (this=0x2c585e0, d=0x224b920, s=1) at /resnick/home/jingze/ATHENAK/builds/inc_cpu_built_in_pgens_dbg/src/src/tasklist/task_list.hpp:93
#23 0x0000000000747230 in TaskList::DoAvailable (this=0x23e1f20, d=0x224b920, s=1) at /resnick/home/jingze/ATHENAK/builds/inc_cpu_built_in_pgens_dbg/src/src/tasklist/task_list.hpp:150
#24 0x000000000073f2ab in Driver::ExecuteTaskList (this=0x224b920, pm=0x25f4410, tl=..., stage=1) at /resnick/home/jingze/ATHENAK/builds/inc_cpu_built_in_pgens_dbg/src/src/driver/driver.cpp:321
#25 0x0000000000741c22 in Driver::Execute (this=0x224b920, pmesh=0x25f4410, pin=0x226f940, pout=0x2c3adf0) at /resnick/home/jingze/ATHENAK/builds/inc_cpu_built_in_pgens_dbg/src/src/driver/driver.cpp:544
#26 0x000000000040e5ee in main (argc=25, argv=0x7fffffffbff8) at /resnick/home/jingze/ATHENAK/builds/inc_cpu_built_in_pgens_dbg/src/src/main.cpp:603
#1  0x00007ffff763e686 in raise () from /lib64/libc.so.6
```

## Binaries

- ATHENA_CPU xthinfix-1009 413b38af, built-in pgens, MPI, Release: `athena_cpu_built_in_pgens_413b38af` md5 `979f0438dab6e2dd29d2dfffc6defdcc`
- ATHENA_GPU 413b38af he_star_m1 H200 nofma (CUDA 12.9, HOPPER90, hpcx, host -ffp-contract=off): `athena_gpu_he_star_m1_413b38af_nofma` md5 `12ce0f6cc19d1f9db258ad74c5c5b5cb`
- ATHENA_BASE rt-integration 5b304cf9 he_star_m1 H200 nofma (reused): `athena_gpu_he_star_m1_5b304cf9_nofma` md5 `a9724eea21f363dc55ef0bc664ab3cf4`

## GPU point (job 4309030, 2 H200; smoke 4308662 sq nlim 3: rc 0, NON-CONVERGED 0, Picard mean 21.3)

Adapted only: #SBATCH header, srun `--mpi=pmix -n 2 -c 8 --cpu-bind=cores gpuwrap.sh`; python with h5py. All 7 arms rc=0 fatal=0.

RESULTS_gpu_arms.txt (raw):
```
base_1    0.463 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.333333e+00', '9.000000e+00', '0.000000e+00')  inner mean 3.330769e+00  plm pos-fallbacks -  FATAL 0
all_1     0.465 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.333333e+00', '9.000000e+00', '0.000000e+00')  inner mean 3.330769e+00  plm pos-fallbacks -  FATAL 0
sq_1      1.901 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.656667e+01', '5.300000e+01', '0.000000e+00')  inner mean 2.947029e+00  plm pos-fallbacks 0.000000e+00  FATAL 0
aq_1      1.825 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.590000e+01', '5.100000e+01', '0.000000e+00')  inner mean 2.946260e+00  plm pos-fallbacks 0.000000e+00  FATAL 0
kn03_1    0.467 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.266667e+00', '8.000000e+00', '0.000000e+00')  inner mean 3.367188e+00  plm pos-fallbacks -  FATAL 0
knp_1     1.752 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.530000e+01', '5.100000e+01', '0.000000e+00')  inner mean 2.902134e+00  plm pos-fallbacks 0.000000e+00  FATAL 0
base_2    0.467 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.333333e+00', '9.000000e+00', '0.000000e+00')  inner mean 3.330769e+00  plm pos-fallbacks -  FATAL 0
-- bin data base_1 vs all_1 (after the parameter header <par_end>), and repeats
base_1 vs all_1: 4 bin files, 0 differ in data
base_1 vs base_2: 4 bin files, 0 differ in data
-- all vs steep / steep+plm by radius (max over angles |a-b|/max|b|), last dump
== all_1 vs sq_1
bands r/R_ph: 0-0.5 0.5-0.8 0.8-0.95 0.95-1.05 1.05-1.3 1.3-2 2-10
agcar3d.hydro_w.00001.bin    dens       0.0e+00 0.0e+00 1.3e-06 1.2e-02 6.0e-03 1.3e-04 5.9e-04
agcar3d.hydro_w.00001.bin    velx       9.8e-08 5.0e-08 5.4e-03 3.9e+00 3.3e-02 2.6e-02 1.0e-01
agcar3d.hydro_w.00001.bin    vely       3.8e-06 5.1e-08 2.8e-04 3.6e-01 9.6e-01 1.1e+00 1.0e+00
agcar3d.hydro_w.00001.bin    velz       5.3e-06 2.2e-08 2.8e-04 3.3e-01 9.7e-01 1.0e+00 1.0e+00
agcar3d.hydro_w.00001.bin    eint       0.0e+00 0.0e+00 4.1e-05 7.6e-02 6.2e-03 8.2e-02 9.5e-02
agcar3d.m1.00001.bin         m1_e       0.0e+00 0.0e+00 1.7e-04 7.6e-02 9.1e-03 2.7e-02 7.5e-02
agcar3d.m1.00001.bin         m1_f1      8.0e-08 2.0e-08 8.5e-05 7.1e-04 1.9e-03 1.3e-02 5.7e-02
agcar3d.m1.00001.bin         m1_f2      5.4e-06 1.8e-08 2.7e-04 2.1e-01 3.8e-01 1.0e+00 1.0e+00
agcar3d.m1.00001.bin         m1_f3      7.9e-06 2.8e-08 2.8e-04 2.0e-01 3.2e-01 9.5e-01 1.1e+00
== all_1 vs aq_1
bands r/R_ph: 0-0.5 0.5-0.8 0.8-0.95 0.95-1.05 1.05-1.3 1.3-2 2-10
agcar3d.hydro_w.00001.bin    dens       0.0e+00 0.0e+00 1.2e-06 1.3e-02 1.5e-03 2.1e-04 1.2e-03
agcar3d.hydro_w.00001.bin    velx       9.8e-08 5.0e-08 5.0e-03 3.1e+00 1.0e-02 3.0e-02 5.1e-02
agcar3d.hydro_w.00001.bin    vely       3.7e-06 3.8e-08 2.8e-04 3.4e-01 1.0e+00 1.0e+00 1.0e+00
agcar3d.hydro_w.00001.bin    velz       3.6e-06 7.0e-09 2.8e-04 3.3e-01 8.7e-01 1.0e+00 1.0e+00
agcar3d.hydro_w.00001.bin    eint       0.0e+00 0.0e+00 3.8e-05 8.2e-02 1.4e-03 8.1e-02 9.5e-02
agcar3d.m1.00001.bin         m1_e       0.0e+00 0.0e+00 1.6e-04 8.1e-02 9.4e-03 2.7e-02 7.4e-02
agcar3d.m1.00001.bin         m1_f1      8.0e-08 3.7e-09 8.0e-05 7.1e-04 2.1e-03 1.4e-02 5.6e-02
agcar3d.m1.00001.bin         m1_f2      5.4e-06 1.8e-08 2.7e-04 2.3e-01 4.9e-01 1.0e+00 1.0e+00
agcar3d.m1.00001.bin         m1_f3      5.4e-06 2.8e-08 2.8e-04 2.2e-01 2.7e-01 1.0e+00 1.1e+00
== all_1 vs kn03_1
bands r/R_ph: 0-0.5 0.5-0.8 0.8-0.95 0.95-1.05 1.05-1.3 1.3-2 2-10
agcar3d.hydro_w.00001.bin    dens       0.0e+00 0.0e+00 2.1e-07 1.3e-03 7.6e-04 2.0e-07 5.5e-07
agcar3d.hydro_w.00001.bin    velx       9.8e-08 3.5e-08 6.3e-04 9.2e-03 5.9e-03 5.4e-05 3.6e-05
agcar3d.hydro_w.00001.bin    vely       3.0e-06 3.8e-08 1.4e-05 2.2e-01 3.3e-01 5.3e-03 5.9e-03
agcar3d.hydro_w.00001.bin    velz       2.9e-06 4.4e-08 1.3e-05 1.6e-01 3.1e-01 8.5e-03 9.4e-03
agcar3d.hydro_w.00001.bin    eint       0.0e+00 0.0e+00 4.9e-06 1.3e-03 7.6e-04 3.1e-05 2.9e-05
agcar3d.m1.00001.bin         m1_e       0.0e+00 0.0e+00 1.9e-05 6.0e-04 2.5e-05 1.5e-05 9.7e-06
agcar3d.m1.00001.bin         m1_f1      7.9e-08 2.0e-08 9.2e-06 1.5e-05 1.3e-05 1.0e-05 8.2e-06
agcar3d.m1.00001.bin         m1_f2      4.4e-06 1.8e-08 1.0e-05 2.1e-02 1.2e-03 6.7e-04 2.3e-04
agcar3d.m1.00001.bin         m1_f3      4.5e-06 2.8e-08 1.1e-05 2.2e-02 7.7e-04 4.8e-04 3.2e-04
== all_1 vs knp_1
bands r/R_ph: 0-0.5 0.5-0.8 0.8-0.95 0.95-1.05 1.05-1.3 1.3-2 2-10
agcar3d.hydro_w.00001.bin    dens       0.0e+00 0.0e+00 1.3e-06 1.3e-02 5.1e-03 1.6e-04 1.1e-03
agcar3d.hydro_w.00001.bin    velx       9.8e-08 5.0e-08 5.6e-03 3.1e+00 2.9e-02 3.0e-02 5.5e-02
agcar3d.hydro_w.00001.bin    vely       3.6e-06 5.1e-08 2.8e-04 3.4e-01 9.5e-01 1.0e+00 1.0e+00
agcar3d.hydro_w.00001.bin    velz       3.3e-06 7.0e-09 2.8e-04 3.1e-01 9.8e-01 1.0e+00 1.0e+00
agcar3d.hydro_w.00001.bin    eint       0.0e+00 0.0e+00 4.3e-05 8.2e-02 5.3e-03 8.3e-02 8.8e-02
agcar3d.m1.00001.bin         m1_e       0.0e+00 0.0e+00 1.7e-04 8.1e-02 9.4e-03 2.7e-02 5.8e-02
agcar3d.m1.00001.bin         m1_f1      8.0e-08 2.0e-08 8.8e-05 6.9e-04 1.8e-03 1.2e-02 4.0e-02
agcar3d.m1.00001.bin         m1_f2      4.9e-06 1.4e-08 2.7e-04 2.0e-01 4.0e-01 1.1e+00 1.1e+00
agcar3d.m1.00001.bin         m1_f3      4.9e-06 2.8e-08 2.7e-04 2.0e-01 2.9e-01 1.0e+00 1.1e+00
```

## RESULTS/order_eval.txt (raw)

```
atm          cen  be       S e      lev 32-64-128-256-512  L1 1.82e-03 4.59e-04 1.18e-04 3.86e-05 | ord  1.99  1.96  1.61
atm          cen  be       S f1     lev 32-64-128-256-512  L1 9.17e-01 7.37e-01 1.12e+00 5.44e+00 | ord  0.32 -0.60 -2.28
atm          hr   be       S e      lev 32-64-128-256-512  L1 2.11e-03 9.83e-04 8.98e-04 1.10e-03 | ord  1.10  0.13 -0.29
atm          hr   be       S f1     lev 32-64-128-256-512  L1 1.27e+00 1.91e+01 5.60e+01 9.75e-01 | ord -3.91 -1.55  5.85
atm          sq   be       S e      lev 32-64-128-256-512  L1 2.06e-03 5.52e-04 1.63e-04 6.14e-05 | ord  1.90  1.76  1.41
atm          sq   be       S f1     lev 32-64-128-256-512  L1 1.47e+00 6.63e+00 1.56e+01 1.34e+00 | ord -2.17 -1.23  3.54
atm          sqf  be       S e      lev 32-64-128-256-512  L1 2.06e-03 5.52e-04 1.63e-04 6.14e-05 | ord  1.90  1.76  1.41
atm          sqf  be       S f1     lev 32-64-128-256-512  L1 1.47e+00 6.63e+00 1.56e+01 1.34e+00 | ord -2.17 -1.23  3.54
atm          sqv  be       S e      lev 32-64-128-256-512  L1 2.06e-03 5.52e-04 1.63e-04 6.14e-05 | ord  1.90  1.76  1.41
atm          sqv  be       S f1     lev 32-64-128-256-512  L1 1.47e+00 6.63e+00 1.56e+01 1.34e+00 | ord -2.17 -1.23  3.54
atm          sqa  be       S e      lev 32-64-128-256-512  L1 2.06e-03 5.52e-04 1.63e-04 6.14e-05 | ord  1.90  1.76  1.41
atm          sqa  be       S f1     lev 32-64-128-256-512  L1 1.47e+00 6.63e+00 1.56e+01 1.34e+00 | ord -2.17 -1.23  3.54
atm          aq   be       S e      lev 32-64-128-256-512  L1 2.08e-03 5.67e-04 1.69e-04 6.41e-05 | ord  1.87  1.74  1.40
atm          aq   be       S f1     lev 32-64-128-256-512  L1 1.51e+00 6.04e+00 1.49e+01 1.29e+00 | ord -2.00 -1.30  3.53
atm          kn01 be       S e      lev 32-64-128-256-512  L1 2.04e-03 8.45e-04 6.26e-04 5.77e-04 | ord  1.27  0.43  0.12
atm          kn01 be       S f1     lev 32-64-128-256-512  L1 1.27e+00 1.88e+01 5.45e+01 1.01e+00 | ord -3.89 -1.54  5.76
atm          kn03 be       S e      lev 32-64-128-256-512  L1 1.95e-03 6.74e-04 3.22e-04 1.33e-04 | ord  1.53  1.07  1.27
atm          kn03 be       S f1     lev 32-64-128-256-512  L1 1.27e+00 1.73e+01 4.54e+01 1.27e+00 | ord -3.77 -1.39  5.16
atm          kn1  be       S e      lev 32-64-128-256-512  L1 1.80e-03 4.70e-04 1.55e-04 1.05e-04 | ord  1.94  1.60  0.56
atm          kn1  be       S f1     lev 32-64-128-256-512  L1 1.46e+00 1.28e+01 1.45e+02 5.74e-01 | ord -3.13 -3.51  7.99
atm          knp  be       S e      lev 32-64-128-256-512  L1 2.08e-03 5.67e-04 1.69e-04 6.40e-05 | ord  1.87  1.74  1.40
atm          knp  be       S f1     lev 32-64-128-256-512  L1 1.51e+00 6.01e+00 1.45e+01 1.33e+00 | ord -1.99 -1.27  3.44
atm          hrx0 be       S e      lev 32-64-128-256-512  L1 1.77e-03 4.41e-04 1.11e-04 3.91e-05 | ord  2.00  1.99  1.51
atm          hrx0 be       S f1     lev 32-64-128-256-512  L1 1.62e+00 1.02e+01 3.29e+01 9.05e-01 | ord -2.65 -1.69  5.18
pulse_k0.128 cen  be       X e      lev 32-64-128-256-512  L1 7.26e-02 2.14e-02 5.49e-03 1.24e-03 | ord  1.76  1.97  2.15
pulse_k0.128 cen  be       X f1     lev 32-64-128-256-512  L1 1.04e-01 2.90e-02 7.32e-03 1.67e-03 | ord  1.83  1.99  2.13
pulse_k0.128 cen  be       T e      lev 1-2-4-8-16  L1 4.31e-02 2.29e-02 1.18e-02 6.02e-03 | ord  0.91  0.95  0.97
pulse_k0.128 cen  be       T f1     lev 1-2-4-8-16  L1 5.40e-02 2.95e-02 1.55e-02 7.95e-03 | ord  0.87  0.93  0.96
pulse_k0.128 cen  be       C e      lev 32-64-128-256-512  L1 8.18e-02 2.68e-02 1.16e-02 5.81e-03 | ord  1.61  1.21  1.00
pulse_k0.128 cen  be       C f1     lev 32-64-128-256-512  L1 1.16e-01 4.63e-02 1.95e-02 8.79e-03 | ord  1.32  1.25  1.15
pulse_k0.128 cen  hesdirk2 X e      lev 32-64-128-256-512  L1 7.49e-02 2.23e-02 5.80e-03 1.41e-03 | ord  1.75  1.94  2.04
pulse_k0.128 cen  hesdirk2 X f1     lev 32-64-128-256-512  L1 1.09e-01 3.06e-02 7.81e-03 1.87e-03 | ord  1.83  1.97  2.06
pulse_k0.128 cen  hesdirk2 T e      lev 1-2-4-8-16  L1 2.38e-03 6.88e-04 2.11e-04 7.09e-05 | ord  1.79  1.71  1.57
pulse_k0.128 cen  hesdirk2 T f1     lev 1-2-4-8-16  L1 3.03e-03 8.44e-04 2.51e-04 8.28e-05 | ord  1.84  1.75  1.60
pulse_k0.128 cen  hesdirk2 C e      lev 32-64-128-256-512  L1 7.62e-02 2.26e-02 5.83e-03 1.41e-03 | ord  1.75  1.95  2.04
pulse_k0.128 cen  hesdirk2 C f1     lev 32-64-128-256-512  L1 1.10e-01 3.09e-02 7.93e-03 1.92e-03 | ord  1.83  1.96  2.05
pulse_k0.128 hr   be       X e      lev 32-64-128-256-512  L1 7.25e-02 2.14e-02 6.33e-03 1.11e-02 | ord  1.76  1.76 -0.82
pulse_k0.128 hr   be       X f1     lev 32-64-128-256-512  L1 1.04e-01 2.91e-02 7.96e-03 1.71e-02 | ord  1.83  1.87 -1.10
pulse_k0.128 hr   be       T e      lev 1-2-4-8-16  L1 5.35e-02 2.81e-02 1.44e-02 7.33e-03 | ord  0.93  0.96  0.98
pulse_k0.128 hr   be       T f1     lev 1-2-4-8-16  L1 6.54e-02 3.51e-02 1.83e-02 9.31e-03 | ord  0.90  0.94  0.97
pulse_k0.128 hr   be       C e      lev 32-64-128-256-512  L1 8.16e-02 2.65e-02 1.07e-02 4.61e-03 | ord  1.62  1.32  1.21
pulse_k0.128 hr   be       C f1     lev 32-64-128-256-512  L1 1.16e-01 4.61e-02 1.90e-02 8.41e-03 | ord  1.33  1.28  1.17
pulse_k0.128 hr   hesdirk2 X e      lev 32-64-128-256-512  L1 7.49e-02 2.23e-02 5.97e-03 3.70e-03 | ord  1.74  1.90  0.69
pulse_k0.128 hr   hesdirk2 X f1     lev 32-64-128-256-512  L1 1.09e-01 3.06e-02 7.88e-03 5.25e-03 | ord  1.83  1.96  0.59
pulse_k0.128 hr   hesdirk2 T e      lev 1-2-4-8-16  L1 5.83e-03 2.47e-03 1.13e-03 5.42e-04 | ord  1.24  1.13  1.06
pulse_k0.128 hr   hesdirk2 T f1     lev 1-2-4-8-16  L1 7.47e-03 3.40e-03 1.64e-03 8.11e-04 | ord  1.14  1.05  1.02
pulse_k0.128 hr   hesdirk2 C e      lev 32-64-128-256-512  L1 7.62e-02 2.26e-02 5.93e-03 1.80e-03 | ord  1.75  1.93  1.72
pulse_k0.128 hr   hesdirk2 C f1     lev 32-64-128-256-512  L1 1.10e-01 3.09e-02 7.95e-03 2.23e-03 | ord  1.83  1.96  1.83
pulse_k0.128 sq   be       X e      lev 32-64-128-256-512  L1 7.26e-02 2.14e-02 5.49e-03 1.24e-03 | ord  1.76  1.97  2.15
pulse_k0.128 sq   be       X f1     lev 32-64-128-256-512  L1 1.04e-01 2.90e-02 7.32e-03 1.67e-03 | ord  1.83  1.99  2.13
pulse_k0.128 sq   be       T e      lev 1-2-4-8-16  L1 4.31e-02 2.29e-02 1.18e-02 6.02e-03 | ord  0.91  0.95  0.97
pulse_k0.128 sq   be       T f1     lev 1-2-4-8-16  L1 5.40e-02 2.95e-02 1.55e-02 7.95e-03 | ord  0.87  0.93  0.96
pulse_k0.128 sq   be       C e      lev 32-64-128-256-512  L1 8.18e-02 2.68e-02 1.16e-02 5.81e-03 | ord  1.61  1.21  1.00
pulse_k0.128 sq   be       C f1     lev 32-64-128-256-512  L1 1.16e-01 4.63e-02 1.95e-02 8.79e-03 | ord  1.32  1.25  1.15
pulse_k0.128 sq   hesdirk2 X e      lev 32-64-128-256-512  L1 7.49e-02 2.23e-02 5.80e-03 1.41e-03 | ord  1.75  1.94  2.04
pulse_k0.128 sq   hesdirk2 X f1     lev 32-64-128-256-512  L1 1.09e-01 3.06e-02 7.81e-03 1.87e-03 | ord  1.83  1.97  2.06
pulse_k0.128 sq   hesdirk2 T e      lev 1-2-4-8-16  L1 2.38e-03 6.88e-04 2.11e-04 7.09e-05 | ord  1.79  1.71  1.57
pulse_k0.128 sq   hesdirk2 T f1     lev 1-2-4-8-16  L1 3.03e-03 8.44e-04 2.51e-04 8.28e-05 | ord  1.84  1.75  1.60
pulse_k0.128 sq   hesdirk2 C e      lev 32-64-128-256-512  L1 7.62e-02 2.26e-02 5.83e-03 1.41e-03 | ord  1.75  1.95  2.04
pulse_k0.128 sq   hesdirk2 C f1     lev 32-64-128-256-512  L1 1.10e-01 3.09e-02 7.93e-03 1.92e-03 | ord  1.83  1.96  2.05
pulse_k0.128 sqf  be       X e      lev 32-64-128-256-512  L1 7.26e-02 2.14e-02 5.49e-03 1.24e-03 | ord  1.76  1.97  2.15
pulse_k0.128 sqf  be       X f1     lev 32-64-128-256-512  L1 1.04e-01 2.90e-02 7.32e-03 1.67e-03 | ord  1.83  1.99  2.13
pulse_k0.128 sqf  be       T e      lev 1-2-4-8-16  L1 4.31e-02 2.29e-02 1.18e-02 6.02e-03 | ord  0.91  0.95  0.97
pulse_k0.128 sqf  be       T f1     lev 1-2-4-8-16  L1 5.40e-02 2.95e-02 1.55e-02 7.95e-03 | ord  0.87  0.93  0.96
pulse_k0.128 sqf  be       C e      lev 32-64-128-256-512  L1 8.18e-02 2.68e-02 1.16e-02 5.81e-03 | ord  1.61  1.21  1.00
pulse_k0.128 sqf  be       C f1     lev 32-64-128-256-512  L1 1.16e-01 4.63e-02 1.95e-02 8.79e-03 | ord  1.32  1.25  1.15
pulse_k0.128 sqf  hesdirk2 X e      lev 32-64-128-256-512  L1 7.49e-02 2.23e-02 5.80e-03 1.41e-03 | ord  1.75  1.94  2.04
pulse_k0.128 sqf  hesdirk2 X f1     lev 32-64-128-256-512  L1 1.09e-01 3.06e-02 7.81e-03 1.87e-03 | ord  1.83  1.97  2.06
pulse_k0.128 sqf  hesdirk2 T e      lev 1-2-4-8-16  L1 2.38e-03 6.88e-04 2.11e-04 7.09e-05 | ord  1.79  1.71  1.57
pulse_k0.128 sqf  hesdirk2 T f1     lev 1-2-4-8-16  L1 3.03e-03 8.44e-04 2.51e-04 8.28e-05 | ord  1.84  1.75  1.60
pulse_k0.128 sqf  hesdirk2 C e      lev 32-64-128-256-512  L1 7.62e-02 2.26e-02 5.83e-03 1.41e-03 | ord  1.75  1.95  2.04
pulse_k0.128 sqf  hesdirk2 C f1     lev 32-64-128-256-512  L1 1.10e-01 3.09e-02 7.93e-03 1.92e-03 | ord  1.83  1.96  2.05
pulse_k0.128 sqv  be       X e      lev 32-64-128-256-512  L1 7.26e-02 2.14e-02 5.49e-03 1.24e-03 | ord  1.76  1.97  2.15
pulse_k0.128 sqv  be       X f1     lev 32-64-128-256-512  L1 1.04e-01 2.90e-02 7.32e-03 1.67e-03 | ord  1.83  1.99  2.13
pulse_k0.128 sqv  be       T e      lev 1-2-4-8-16  L1 4.31e-02 2.29e-02 1.18e-02 6.02e-03 | ord  0.91  0.95  0.97
pulse_k0.128 sqv  be       T f1     lev 1-2-4-8-16  L1 5.40e-02 2.95e-02 1.55e-02 7.95e-03 | ord  0.87  0.93  0.96
pulse_k0.128 sqv  be       C e      lev 32-64-128-256-512  L1 8.18e-02 2.68e-02 1.16e-02 5.81e-03 | ord  1.61  1.21  1.00
pulse_k0.128 sqv  be       C f1     lev 32-64-128-256-512  L1 1.16e-01 4.63e-02 1.95e-02 8.79e-03 | ord  1.32  1.25  1.15
pulse_k0.128 sqv  hesdirk2 X e      lev 32-64-128-256-512  L1 7.49e-02 2.23e-02 5.80e-03 1.41e-03 | ord  1.75  1.94  2.04
pulse_k0.128 sqv  hesdirk2 X f1     lev 32-64-128-256-512  L1 1.09e-01 3.06e-02 7.81e-03 1.87e-03 | ord  1.83  1.97  2.06
pulse_k0.128 sqv  hesdirk2 T e      lev 1-2-4-8-16  L1 2.38e-03 6.88e-04 2.11e-04 7.09e-05 | ord  1.79  1.71  1.57
pulse_k0.128 sqv  hesdirk2 T f1     lev 1-2-4-8-16  L1 3.03e-03 8.44e-04 2.51e-04 8.28e-05 | ord  1.84  1.75  1.60
pulse_k0.128 sqv  hesdirk2 C e      lev 32-64-128-256-512  L1 7.62e-02 2.26e-02 5.83e-03 1.41e-03 | ord  1.75  1.95  2.04
pulse_k0.128 sqv  hesdirk2 C f1     lev 32-64-128-256-512  L1 1.10e-01 3.09e-02 7.93e-03 1.92e-03 | ord  1.83  1.96  2.05
pulse_k0.128 sqa  be       X e      lev 32-64-128-256-512  L1 7.26e-02 2.14e-02 5.49e-03 1.24e-03 | ord  1.76  1.97  2.15
pulse_k0.128 sqa  be       X f1     lev 32-64-128-256-512  L1 1.04e-01 2.90e-02 7.32e-03 1.67e-03 | ord  1.83  1.99  2.13
pulse_k0.128 sqa  be       T e      lev 1-2-4-8-16  L1 4.31e-02 2.29e-02 1.18e-02 6.02e-03 | ord  0.91  0.95  0.97
pulse_k0.128 sqa  be       T f1     lev 1-2-4-8-16  L1 5.40e-02 2.95e-02 1.55e-02 7.95e-03 | ord  0.87  0.93  0.96
pulse_k0.128 sqa  be       C e      lev 32-64-128-256-512  L1 8.18e-02 2.68e-02 1.16e-02 5.81e-03 | ord  1.61  1.21  1.00
pulse_k0.128 sqa  be       C f1     lev 32-64-128-256-512  L1 1.16e-01 4.63e-02 1.95e-02 8.79e-03 | ord  1.32  1.25  1.15
pulse_k0.128 sqa  hesdirk2 X e      lev 32-64-128-256-512  L1 7.49e-02 2.23e-02 5.80e-03 1.41e-03 | ord  1.75  1.94  2.04
pulse_k0.128 sqa  hesdirk2 X f1     lev 32-64-128-256-512  L1 1.09e-01 3.06e-02 7.81e-03 1.87e-03 | ord  1.83  1.97  2.06
pulse_k0.128 sqa  hesdirk2 T e      lev 1-2-4-8-16  L1 2.38e-03 6.88e-04 2.11e-04 7.09e-05 | ord  1.79  1.71  1.57
pulse_k0.128 sqa  hesdirk2 T f1     lev 1-2-4-8-16  L1 3.03e-03 8.44e-04 2.51e-04 8.28e-05 | ord  1.84  1.75  1.60
pulse_k0.128 sqa  hesdirk2 C e      lev 32-64-128-256-512  L1 7.62e-02 2.26e-02 5.83e-03 1.41e-03 | ord  1.75  1.95  2.04
pulse_k0.128 sqa  hesdirk2 C f1     lev 32-64-128-256-512  L1 1.10e-01 3.09e-02 7.93e-03 1.92e-03 | ord  1.83  1.96  2.05
pulse_k0.128 aq   be       X e      lev 32-64-128-256-512  L1 7.25e-02 2.14e-02 6.30e-03 1.12e-02 | ord  1.76  1.76 -0.83
pulse_k0.128 aq   be       X f1     lev 32-64-128-256-512  L1 1.04e-01 2.90e-02 7.92e-03 1.71e-02 | ord  1.83  1.87 -1.11
pulse_k0.128 aq   be       T e      lev 1-2-4-8-16  L1 5.35e-02 2.81e-02 1.45e-02 7.34e-03 | ord  0.93  0.96  0.98
pulse_k0.128 aq   be       T f1     lev 1-2-4-8-16  L1 6.55e-02 3.52e-02 1.83e-02 9.32e-03 | ord  0.90  0.94  0.97
pulse_k0.128 aq   be       C e      lev 32-64-128-256-512  L1 8.15e-02 2.65e-02 1.06e-02 4.61e-03 | ord  1.62  1.32  1.21
pulse_k0.128 aq   be       C f1     lev 32-64-128-256-512  L1 1.16e-01 4.61e-02 1.90e-02 8.40e-03 | ord  1.33  1.28  1.17
pulse_k0.128 aq   hesdirk2 X e      lev 32-64-128-256-512  L1 7.49e-02 2.23e-02 5.96e-03 3.72e-03 | ord  1.74  1.91  0.68
pulse_k0.128 aq   hesdirk2 X f1     lev 32-64-128-256-512  L1 1.09e-01 3.06e-02 7.87e-03 5.28e-03 | ord  1.83  1.96  0.58
pulse_k0.128 aq   hesdirk2 T e      lev 1-2-4-8-16  L1 5.92e-03 2.51e-03 1.16e-03 5.54e-04 | ord  1.23  1.12  1.06
pulse_k0.128 aq   hesdirk2 T f1     lev 1-2-4-8-16  L1 7.52e-03 3.43e-03 1.66e-03 8.20e-04 | ord  1.13  1.04  1.02
pulse_k0.128 aq   hesdirk2 C e      lev 32-64-128-256-512  L1 7.61e-02 2.26e-02 5.93e-03 1.80e-03 | ord  1.75  1.93  1.72
pulse_k0.128 aq   hesdirk2 C f1     lev 32-64-128-256-512  L1 1.10e-01 3.09e-02 7.95e-03 2.23e-03 | ord  1.83  1.96  1.83
pulse_k0.128 kn01 be       X e      lev 32-64-128-256-512  L1 7.25e-02 2.14e-02 6.31e-03 1.12e-02 | ord  1.76  1.76 -0.82
pulse_k0.128 kn01 be       X f1     lev 32-64-128-256-512  L1 1.04e-01 2.90e-02 7.52e-03 1.50e-02 | ord  1.83  1.95 -1.00
pulse_k0.128 kn01 be       T e      lev 1-2-4-8-16  L1 5.36e-02 2.82e-02 1.44e-02 7.33e-03 | ord  0.93  0.96  0.98
pulse_k0.128 kn01 be       T f1     lev 1-2-4-8-16  L1 6.33e-02 3.40e-02 1.76e-02 9.00e-03 | ord  0.90  0.95  0.97
pulse_k0.128 kn01 be       C e      lev 32-64-128-256-512  L1 8.16e-02 2.65e-02 1.06e-02 4.59e-03 | ord  1.62  1.32  1.21
pulse_k0.128 kn01 be       C f1     lev 32-64-128-256-512  L1 1.16e-01 4.61e-02 1.87e-02 7.79e-03 | ord  1.33  1.30  1.27
pulse_k0.128 kn01 hesdirk2 X e      lev 32-64-128-256-512  L1 7.49e-02 2.23e-02 5.97e-03 3.69e-03 | ord  1.74  1.90  0.69
pulse_k0.128 kn01 hesdirk2 X f1     lev 32-64-128-256-512  L1 1.09e-01 3.06e-02 7.81e-03 4.80e-03 | ord  1.83  1.97  0.70
pulse_k0.128 kn01 hesdirk2 T e      lev 1-2-4-8-16  L1 5.86e-03 2.49e-03 1.14e-03 5.49e-04 | ord  1.24  1.12  1.06
pulse_k0.128 kn01 hesdirk2 T f1     lev 1-2-4-8-16  L1 6.72e-03 3.05e-03 1.47e-03 7.27e-04 | ord  1.14  1.05  1.02
pulse_k0.128 kn01 hesdirk2 C e      lev 32-64-128-256-512  L1 7.62e-02 2.26e-02 5.92e-03 1.80e-03 | ord  1.75  1.93  1.72
pulse_k0.128 kn01 hesdirk2 C f1     lev 32-64-128-256-512  L1 1.10e-01 3.09e-02 7.92e-03 2.11e-03 | ord  1.83  1.96  1.91
pulse_k0.128 kn03 be       X e      lev 32-64-128-256-512  L1 7.25e-02 2.14e-02 6.30e-03 1.11e-02 | ord  1.76  1.77 -0.81
pulse_k0.128 kn03 be       X f1     lev 32-64-128-256-512  L1 1.04e-01 2.90e-02 7.52e-03 1.49e-02 | ord  1.83  1.95 -0.99
pulse_k0.128 kn03 be       T e      lev 1-2-4-8-16  L1 5.35e-02 2.81e-02 1.44e-02 7.32e-03 | ord  0.93  0.96  0.98
pulse_k0.128 kn03 be       T f1     lev 1-2-4-8-16  L1 6.29e-02 3.38e-02 1.76e-02 8.98e-03 | ord  0.89  0.94  0.97
pulse_k0.128 kn03 be       C e      lev 32-64-128-256-512  L1 8.16e-02 2.65e-02 1.06e-02 4.58e-03 | ord  1.62  1.32  1.22
pulse_k0.128 kn03 be       C f1     lev 32-64-128-256-512  L1 1.16e-01 4.61e-02 1.87e-02 7.78e-03 | ord  1.33  1.30  1.27
pulse_k0.128 kn03 hesdirk2 X e      lev 32-64-128-256-512  L1 7.49e-02 2.23e-02 5.97e-03 3.67e-03 | ord  1.74  1.90  0.70
pulse_k0.128 kn03 hesdirk2 X f1     lev 32-64-128-256-512  L1 1.09e-01 3.06e-02 7.81e-03 4.77e-03 | ord  1.83  1.97  0.71
pulse_k0.128 kn03 hesdirk2 T e      lev 1-2-4-8-16  L1 5.84e-03 2.47e-03 1.13e-03 5.45e-04 | ord  1.24  1.12  1.06
pulse_k0.128 kn03 hesdirk2 T f1     lev 1-2-4-8-16  L1 6.63e-03 3.00e-03 1.45e-03 7.18e-04 | ord  1.14  1.05  1.02
pulse_k0.128 kn03 hesdirk2 C e      lev 32-64-128-256-512  L1 7.62e-02 2.26e-02 5.92e-03 1.80e-03 | ord  1.75  1.93  1.72
pulse_k0.128 kn03 hesdirk2 C f1     lev 32-64-128-256-512  L1 1.10e-01 3.09e-02 7.92e-03 2.11e-03 | ord  1.83  1.96  1.91
pulse_k0.128 kn1  be       X e      lev 32-64-128-256-512  L1 7.25e-02 2.14e-02 6.29e-03 1.09e-02 | ord  1.76  1.77 -0.79
pulse_k0.128 kn1  be       X f1     lev 32-64-128-256-512  L1 1.04e-01 2.90e-02 7.51e-03 1.47e-02 | ord  1.84  1.95 -0.97
pulse_k0.128 kn1  be       T e      lev 1-2-4-8-16  L1 5.32e-02 2.80e-02 1.44e-02 7.31e-03 | ord  0.92  0.96  0.98
pulse_k0.128 kn1  be       T f1     lev 1-2-4-8-16  L1 6.24e-02 3.37e-02 1.75e-02 8.96e-03 | ord  0.89  0.94  0.97
pulse_k0.128 kn1  be       C e      lev 32-64-128-256-512  L1 8.16e-02 2.65e-02 1.07e-02 4.59e-03 | ord  1.62  1.32  1.22
pulse_k0.128 kn1  be       C f1     lev 32-64-128-256-512  L1 1.16e-01 4.61e-02 1.88e-02 7.81e-03 | ord  1.33  1.30  1.26
pulse_k0.128 kn1  hesdirk2 X e      lev 32-64-128-256-512  L1 7.49e-02 2.23e-02 5.96e-03 3.63e-03 | ord  1.74  1.91  0.72
pulse_k0.128 kn1  hesdirk2 X f1     lev 32-64-128-256-512  L1 1.09e-01 3.06e-02 7.81e-03 4.68e-03 | ord  1.83  1.97  0.74
pulse_k0.128 kn1  hesdirk2 T e      lev 1-2-4-8-16  L1 5.78e-03 2.44e-03 1.12e-03 5.37e-04 | ord  1.24  1.13  1.06
pulse_k0.128 kn1  hesdirk2 T f1     lev 1-2-4-8-16  L1 6.54e-03 2.94e-03 1.41e-03 6.97e-04 | ord  1.15  1.06  1.02
pulse_k0.128 kn1  hesdirk2 C e      lev 32-64-128-256-512  L1 7.62e-02 2.26e-02 5.92e-03 1.79e-03 | ord  1.75  1.93  1.73
pulse_k0.128 kn1  hesdirk2 C f1     lev 32-64-128-256-512  L1 1.10e-01 3.09e-02 7.92e-03 2.10e-03 | ord  1.83  1.96  1.92
pulse_k0.128 knp  be       X e      lev 32-64-128-256-512  L1 7.25e-02 2.14e-02 6.28e-03 1.11e-02 | ord  1.76  1.77 -0.83
pulse_k0.128 knp  be       X f1     lev 32-64-128-256-512  L1 1.04e-01 2.90e-02 7.48e-03 1.50e-02 | ord  1.84  1.95 -1.00
pulse_k0.128 knp  be       T e      lev 1-2-4-8-16  L1 5.35e-02 2.81e-02 1.44e-02 7.33e-03 | ord  0.93  0.96  0.98
pulse_k0.128 knp  be       T f1     lev 1-2-4-8-16  L1 6.30e-02 3.39e-02 1.76e-02 9.00e-03 | ord  0.89  0.94  0.97
pulse_k0.128 knp  be       C e      lev 32-64-128-256-512  L1 8.15e-02 2.65e-02 1.06e-02 4.58e-03 | ord  1.62  1.32  1.22
pulse_k0.128 knp  be       C f1     lev 32-64-128-256-512  L1 1.16e-01 4.60e-02 1.87e-02 7.78e-03 | ord  1.33  1.30  1.27
pulse_k0.128 knp  hesdirk2 X e      lev 32-64-128-256-512  L1 7.49e-02 2.23e-02 5.95e-03 3.69e-03 | ord  1.74  1.91  0.69
pulse_k0.128 knp  hesdirk2 X f1     lev 32-64-128-256-512  L1 1.09e-01 3.06e-02 7.80e-03 4.79e-03 | ord  1.83  1.97  0.70
pulse_k0.128 knp  hesdirk2 T e      lev 1-2-4-8-16  L1 5.93e-03 2.52e-03 1.16e-03 5.57e-04 | ord  1.24  1.12  1.06
pulse_k0.128 knp  hesdirk2 T f1     lev 1-2-4-8-16  L1 6.68e-03 3.03e-03 1.47e-03 7.27e-04 | ord  1.14  1.04  1.02
pulse_k0.128 knp  hesdirk2 C e      lev 32-64-128-256-512  L1 7.62e-02 2.26e-02 5.92e-03 1.80e-03 | ord  1.75  1.93  1.72
pulse_k0.128 knp  hesdirk2 C f1     lev 32-64-128-256-512  L1 1.10e-01 3.09e-02 7.92e-03 2.11e-03 | ord  1.83  1.96  1.91
pulse_k12.8  cen  be       X e      lev 32-64-128-256-512  L1 2.36e-02 6.37e-03 1.33e-03 2.80e-04 | ord  1.89  2.25  2.25
pulse_k12.8  cen  be       X f1     lev 32-64-128-256-512  L1 5.56e-02 1.63e-02 3.58e-03 7.68e-04 | ord  1.77  2.18  2.22
pulse_k12.8  cen  be       T e      lev 1-2-4-8-16  L1 1.37e-02 7.28e-03 3.77e-03 1.92e-03 | ord  0.91  0.95  0.97
pulse_k12.8  cen  be       T f1     lev 1-2-4-8-16  L1 3.21e-02 1.77e-02 9.29e-03 4.77e-03 | ord  0.86  0.93  0.96
pulse_k12.8  cen  be       C e      lev 32-64-128-256-512  L1 2.81e-02 1.17e-02 4.51e-03 2.04e-03 | ord  1.27  1.37  1.15
pulse_k12.8  cen  be       C f1     lev 32-64-128-256-512  L1 7.55e-02 3.18e-02 1.25e-02 5.48e-03 | ord  1.24  1.35  1.19
pulse_k12.8  cen  hesdirk2 X e      lev 32-64-128-256-512  L1 2.44e-02 6.51e-03 1.37e-03 2.89e-04 | ord  1.90  2.25  2.24
pulse_k12.8  cen  hesdirk2 X f1     lev 32-64-128-256-512  L1 5.69e-02 1.68e-02 3.71e-03 7.99e-04 | ord  1.76  2.18  2.22
pulse_k12.8  cen  hesdirk2 T e      lev 1-2-4-8-16  L1 3.83e-03 1.85e-03 9.03e-04 4.47e-04 | ord  1.05  1.03  1.02
pulse_k12.8  cen  hesdirk2 T f1     lev 1-2-4-8-16  L1 8.30e-03 4.06e-03 2.00e-03 9.93e-04 | ord  1.03  1.02  1.01
pulse_k12.8  cen  hesdirk2 C e      lev 32-64-128-256-512  L1 2.36e-02 7.37e-03 1.91e-03 6.01e-04 | ord  1.68  1.95  1.67
pulse_k12.8  cen  hesdirk2 C f1     lev 32-64-128-256-512  L1 6.12e-02 1.98e-02 5.35e-03 1.69e-03 | ord  1.63  1.89  1.67
pulse_k12.8  hr   be       X e      lev 32-64-128-256-512  L1 2.36e-02 6.36e-03 1.41e-03 2.05e-03 | ord  1.89  2.18 -0.55
pulse_k12.8  hr   be       X f1     lev 32-64-128-256-512  L1 5.56e-02 1.62e-02 3.49e-03 4.78e-03 | ord  1.78  2.21 -0.45
pulse_k12.8  hr   be       T e      lev 1-2-4-8-16  L1 1.54e-02 8.30e-03 4.32e-03 2.21e-03 | ord  0.89  0.94  0.97
pulse_k12.8  hr   be       T f1     lev 1-2-4-8-16  L1 3.43e-02 1.91e-02 1.01e-02 5.19e-03 | ord  0.85  0.92  0.96
pulse_k12.8  hr   be       C e      lev 32-64-128-256-512  L1 2.80e-02 1.15e-02 4.24e-03 1.56e-03 | ord  1.28  1.44  1.44
pulse_k12.8  hr   be       C f1     lev 32-64-128-256-512  L1 7.53e-02 3.17e-02 1.23e-02 5.06e-03 | ord  1.25  1.37  1.28
pulse_k12.8  hr   hesdirk2 X e      lev 32-64-128-256-512  L1 2.44e-02 6.51e-03 1.38e-03 6.57e-04 | ord  1.90  2.24  1.07
pulse_k12.8  hr   hesdirk2 X f1     lev 32-64-128-256-512  L1 5.69e-02 1.67e-02 3.65e-03 1.37e-03 | ord  1.76  2.20  1.41
pulse_k12.8  hr   hesdirk2 T e      lev 1-2-4-8-16  L1 4.39e-03 2.16e-03 1.07e-03 5.30e-04 | ord  1.03  1.02  1.01
pulse_k12.8  hr   hesdirk2 T f1     lev 1-2-4-8-16  L1 9.26e-03 4.58e-03 2.27e-03 1.13e-03 | ord  1.02  1.01  1.01
pulse_k12.8  hr   hesdirk2 C e      lev 32-64-128-256-512  L1 2.36e-02 7.34e-03 1.86e-03 4.83e-04 | ord  1.69  1.98  1.94
pulse_k12.8  hr   hesdirk2 C f1     lev 32-64-128-256-512  L1 6.12e-02 1.98e-02 5.26e-03 1.56e-03 | ord  1.63  1.91  1.76
pulse_k12.8  sq   be       X e      lev 32-64-128-256-512  L1 2.36e-02 6.37e-03 1.33e-03 2.80e-04 | ord  1.89  2.25  2.25
pulse_k12.8  sq   be       X f1     lev 32-64-128-256-512  L1 5.56e-02 1.63e-02 3.58e-03 7.68e-04 | ord  1.77  2.18  2.22
pulse_k12.8  sq   be       T e      lev 1-2-4-8-16  L1 1.37e-02 7.28e-03 3.77e-03 1.92e-03 | ord  0.91  0.95  0.97
pulse_k12.8  sq   be       T f1     lev 1-2-4-8-16  L1 3.21e-02 1.77e-02 9.29e-03 4.77e-03 | ord  0.86  0.93  0.96
pulse_k12.8  sq   be       C e      lev 32-64-128-256-512  L1 2.81e-02 1.17e-02 4.51e-03 2.04e-03 | ord  1.27  1.37  1.15
pulse_k12.8  sq   be       C f1     lev 32-64-128-256-512  L1 7.55e-02 3.18e-02 1.25e-02 5.48e-03 | ord  1.24  1.35  1.19
pulse_k12.8  sq   hesdirk2 X e      lev 32-64-128-256-512  L1 2.44e-02 6.51e-03 1.37e-03 2.89e-04 | ord  1.90  2.25  2.24
pulse_k12.8  sq   hesdirk2 X f1     lev 32-64-128-256-512  L1 5.69e-02 1.68e-02 3.71e-03 7.99e-04 | ord  1.76  2.18  2.22
pulse_k12.8  sq   hesdirk2 T e      lev 1-2-4-8-16  L1 3.83e-03 1.85e-03 9.03e-04 4.47e-04 | ord  1.05  1.03  1.02
pulse_k12.8  sq   hesdirk2 T f1     lev 1-2-4-8-16  L1 8.30e-03 4.06e-03 2.00e-03 9.93e-04 | ord  1.03  1.02  1.01
pulse_k12.8  sq   hesdirk2 C e      lev 32-64-128-256-512  L1 2.36e-02 7.37e-03 1.91e-03 6.01e-04 | ord  1.68  1.95  1.67
pulse_k12.8  sq   hesdirk2 C f1     lev 32-64-128-256-512  L1 6.12e-02 1.98e-02 5.35e-03 1.69e-03 | ord  1.63  1.89  1.67
pulse_k12.8  sqf  be       X e      lev 32-64-128-256-512  L1 2.36e-02 6.37e-03 1.33e-03 2.80e-04 | ord  1.89  2.25  2.25
pulse_k12.8  sqf  be       X f1     lev 32-64-128-256-512  L1 5.56e-02 1.63e-02 3.58e-03 7.68e-04 | ord  1.77  2.18  2.22
pulse_k12.8  sqf  be       T e      lev 1-2-4-8-16  L1 1.37e-02 7.28e-03 3.77e-03 1.92e-03 | ord  0.91  0.95  0.97
pulse_k12.8  sqf  be       T f1     lev 1-2-4-8-16  L1 3.21e-02 1.77e-02 9.29e-03 4.77e-03 | ord  0.86  0.93  0.96
pulse_k12.8  sqf  be       C e      lev 32-64-128-256-512  L1 2.81e-02 1.17e-02 4.51e-03 2.04e-03 | ord  1.27  1.37  1.15
pulse_k12.8  sqf  be       C f1     lev 32-64-128-256-512  L1 7.55e-02 3.18e-02 1.25e-02 5.48e-03 | ord  1.24  1.35  1.19
pulse_k12.8  sqf  hesdirk2 X e      lev 32-64-128-256-512  L1 2.44e-02 6.51e-03 1.37e-03 2.89e-04 | ord  1.90  2.25  2.24
pulse_k12.8  sqf  hesdirk2 X f1     lev 32-64-128-256-512  L1 5.69e-02 1.68e-02 3.71e-03 7.99e-04 | ord  1.76  2.18  2.22
pulse_k12.8  sqf  hesdirk2 T e      lev 1-2-4-8-16  L1 3.83e-03 1.85e-03 9.03e-04 4.47e-04 | ord  1.05  1.03  1.02
pulse_k12.8  sqf  hesdirk2 T f1     lev 1-2-4-8-16  L1 8.30e-03 4.06e-03 2.00e-03 9.93e-04 | ord  1.03  1.02  1.01
pulse_k12.8  sqf  hesdirk2 C e      lev 32-64-128-256-512  L1 2.36e-02 7.37e-03 1.91e-03 6.01e-04 | ord  1.68  1.95  1.67
pulse_k12.8  sqf  hesdirk2 C f1     lev 32-64-128-256-512  L1 6.12e-02 1.98e-02 5.35e-03 1.69e-03 | ord  1.63  1.89  1.67
pulse_k12.8  sqv  be       X e      lev 32-64-128-256-512  L1 2.36e-02 6.37e-03 1.33e-03 2.80e-04 | ord  1.89  2.25  2.25
pulse_k12.8  sqv  be       X f1     lev 32-64-128-256-512  L1 5.56e-02 1.63e-02 3.58e-03 7.68e-04 | ord  1.77  2.18  2.22
pulse_k12.8  sqv  be       T e      lev 1-2-4-8-16  L1 1.37e-02 7.28e-03 3.77e-03 1.92e-03 | ord  0.91  0.95  0.97
pulse_k12.8  sqv  be       T f1     lev 1-2-4-8-16  L1 3.21e-02 1.77e-02 9.29e-03 4.77e-03 | ord  0.86  0.93  0.96
pulse_k12.8  sqv  be       C e      lev 32-64-128-256-512  L1 2.81e-02 1.17e-02 4.51e-03 2.04e-03 | ord  1.27  1.37  1.15
pulse_k12.8  sqv  be       C f1     lev 32-64-128-256-512  L1 7.55e-02 3.18e-02 1.25e-02 5.48e-03 | ord  1.24  1.35  1.19
pulse_k12.8  sqv  hesdirk2 X e      lev 32-64-128-256-512  L1 2.44e-02 6.51e-03 1.37e-03 2.89e-04 | ord  1.90  2.25  2.24
pulse_k12.8  sqv  hesdirk2 X f1     lev 32-64-128-256-512  L1 5.69e-02 1.68e-02 3.71e-03 7.99e-04 | ord  1.76  2.18  2.22
pulse_k12.8  sqv  hesdirk2 T e      lev 1-2-4-8-16  L1 3.83e-03 1.85e-03 9.03e-04 4.47e-04 | ord  1.05  1.03  1.02
pulse_k12.8  sqv  hesdirk2 T f1     lev 1-2-4-8-16  L1 8.30e-03 4.06e-03 2.00e-03 9.93e-04 | ord  1.03  1.02  1.01
pulse_k12.8  sqv  hesdirk2 C e      lev 32-64-128-256-512  L1 2.36e-02 7.37e-03 1.91e-03 6.01e-04 | ord  1.68  1.95  1.67
pulse_k12.8  sqv  hesdirk2 C f1     lev 32-64-128-256-512  L1 6.12e-02 1.98e-02 5.35e-03 1.69e-03 | ord  1.63  1.89  1.67
pulse_k12.8  sqa  be       X e      lev 32-64-128-256-512  L1 2.36e-02 6.37e-03 1.33e-03 2.80e-04 | ord  1.89  2.25  2.25
pulse_k12.8  sqa  be       X f1     lev 32-64-128-256-512  L1 5.56e-02 1.63e-02 3.58e-03 7.68e-04 | ord  1.77  2.18  2.22
pulse_k12.8  sqa  be       T e      lev 1-2-4-8-16  L1 1.37e-02 7.28e-03 3.77e-03 1.92e-03 | ord  0.91  0.95  0.97
pulse_k12.8  sqa  be       T f1     lev 1-2-4-8-16  L1 3.21e-02 1.77e-02 9.29e-03 4.77e-03 | ord  0.86  0.93  0.96
pulse_k12.8  sqa  be       C e      lev 32-64-128-256-512  L1 2.81e-02 1.17e-02 4.51e-03 2.04e-03 | ord  1.27  1.37  1.15
pulse_k12.8  sqa  be       C f1     lev 32-64-128-256-512  L1 7.55e-02 3.18e-02 1.25e-02 5.48e-03 | ord  1.24  1.35  1.19
pulse_k12.8  sqa  hesdirk2 X e      lev 32-64-128-256-512  L1 2.44e-02 6.51e-03 1.37e-03 2.89e-04 | ord  1.90  2.25  2.24
pulse_k12.8  sqa  hesdirk2 X f1     lev 32-64-128-256-512  L1 5.69e-02 1.68e-02 3.71e-03 7.99e-04 | ord  1.76  2.18  2.22
pulse_k12.8  sqa  hesdirk2 T e      lev 1-2-4-8-16  L1 3.83e-03 1.85e-03 9.03e-04 4.47e-04 | ord  1.05  1.03  1.02
pulse_k12.8  sqa  hesdirk2 T f1     lev 1-2-4-8-16  L1 8.30e-03 4.06e-03 2.00e-03 9.93e-04 | ord  1.03  1.02  1.01
pulse_k12.8  sqa  hesdirk2 C e      lev 32-64-128-256-512  L1 2.36e-02 7.37e-03 1.91e-03 6.01e-04 | ord  1.68  1.95  1.67
pulse_k12.8  sqa  hesdirk2 C f1     lev 32-64-128-256-512  L1 6.12e-02 1.98e-02 5.35e-03 1.69e-03 | ord  1.63  1.89  1.67
pulse_k12.8  aq   be       X e      lev 32-64-128-256-512  L1 2.36e-02 6.34e-03 1.39e-03 2.08e-03 | ord  1.89  2.19 -0.58
pulse_k12.8  aq   be       X f1     lev 32-64-128-256-512  L1 5.56e-02 1.61e-02 3.43e-03 4.77e-03 | ord  1.78  2.23 -0.47
pulse_k12.8  aq   be       T e      lev 1-2-4-8-16  L1 1.55e-02 8.37e-03 4.36e-03 2.23e-03 | ord  0.89  0.94  0.97
pulse_k12.8  aq   be       T f1     lev 1-2-4-8-16  L1 3.46e-02 1.92e-02 1.02e-02 5.23e-03 | ord  0.85  0.92  0.96
pulse_k12.8  aq   be       C e      lev 32-64-128-256-512  L1 2.80e-02 1.15e-02 4.24e-03 1.56e-03 | ord  1.28  1.44  1.44
pulse_k12.8  aq   be       C f1     lev 32-64-128-256-512  L1 7.53e-02 3.17e-02 1.22e-02 5.06e-03 | ord  1.25  1.37  1.28
pulse_k12.8  aq   hesdirk2 X e      lev 32-64-128-256-512  L1 2.44e-02 6.50e-03 1.37e-03 6.56e-04 | ord  1.90  2.24  1.07
pulse_k12.8  aq   hesdirk2 X f1     lev 32-64-128-256-512  L1 5.68e-02 1.67e-02 3.63e-03 1.35e-03 | ord  1.76  2.21  1.43
pulse_k12.8  aq   hesdirk2 T e      lev 1-2-4-8-16  L1 4.44e-03 2.18e-03 1.08e-03 5.35e-04 | ord  1.03  1.02  1.01
pulse_k12.8  aq   hesdirk2 T f1     lev 1-2-4-8-16  L1 9.32e-03 4.61e-03 2.29e-03 1.14e-03 | ord  1.02  1.01  1.01
pulse_k12.8  aq   hesdirk2 C e      lev 32-64-128-256-512  L1 2.36e-02 7.34e-03 1.86e-03 4.83e-04 | ord  1.69  1.98  1.94
pulse_k12.8  aq   hesdirk2 C f1     lev 32-64-128-256-512  L1 6.12e-02 1.98e-02 5.26e-03 1.56e-03 | ord  1.63  1.91  1.76
pulse_k12.8  kn01 be       X e      lev 32-64-128-256-512  L1 2.36e-02 6.35e-03 1.29e-03 1.19e-03 | ord  1.89  2.30  0.12
pulse_k12.8  kn01 be       X f1     lev 32-64-128-256-512  L1 5.56e-02 1.62e-02 3.22e-03 2.95e-03 | ord  1.78  2.33  0.13
pulse_k12.8  kn01 be       T e      lev 1-2-4-8-16  L1 1.48e-02 7.92e-03 4.11e-03 2.10e-03 | ord  0.90  0.95  0.97
pulse_k12.8  kn01 be       T f1     lev 1-2-4-8-16  L1 3.42e-02 1.89e-02 1.00e-02 5.14e-03 | ord  0.85  0.92  0.96
pulse_k12.8  kn01 be       C e      lev 32-64-128-256-512  L1 2.80e-02 1.16e-02 4.32e-03 1.64e-03 | ord  1.28  1.42  1.40
pulse_k12.8  kn01 be       C f1     lev 32-64-128-256-512  L1 7.54e-02 3.17e-02 1.21e-02 4.69e-03 | ord  1.25  1.38  1.37
pulse_k12.8  kn01 hesdirk2 X e      lev 32-64-128-256-512  L1 2.44e-02 6.51e-03 1.35e-03 4.29e-04 | ord  1.90  2.27  1.65
pulse_k12.8  kn01 hesdirk2 X f1     lev 32-64-128-256-512  L1 5.69e-02 1.67e-02 3.60e-03 8.93e-04 | ord  1.76  2.22  2.01
pulse_k12.8  kn01 hesdirk2 T e      lev 1-2-4-8-16  L1 4.16e-03 2.03e-03 1.00e-03 4.98e-04 | ord  1.03  1.02  1.01
pulse_k12.8  kn01 hesdirk2 T f1     lev 1-2-4-8-16  L1 9.13e-03 4.50e-03 2.23e-03 1.11e-03 | ord  1.02  1.01  1.01
pulse_k12.8  kn01 hesdirk2 C e      lev 32-64-128-256-512  L1 2.36e-02 7.34e-03 1.86e-03 4.88e-04 | ord  1.68  1.98  1.93
pulse_k12.8  kn01 hesdirk2 C f1     lev 32-64-128-256-512  L1 6.12e-02 1.98e-02 5.23e-03 1.45e-03 | ord  1.63  1.92  1.85
pulse_k12.8  kn03 be       X e      lev 32-64-128-256-512  L1 2.36e-02 6.36e-03 1.28e-03 6.70e-04 | ord  1.89  2.31  0.94
pulse_k12.8  kn03 be       X f1     lev 32-64-128-256-512  L1 5.56e-02 1.62e-02 3.33e-03 1.86e-03 | ord  1.78  2.28  0.84
pulse_k12.8  kn03 be       T e      lev 1-2-4-8-16  L1 1.41e-02 7.58e-03 3.94e-03 2.01e-03 | ord  0.90  0.94  0.97
pulse_k12.8  kn03 be       T f1     lev 1-2-4-8-16  L1 3.31e-02 1.83e-02 9.67e-03 4.97e-03 | ord  0.85  0.92  0.96
pulse_k12.8  kn03 be       C e      lev 32-64-128-256-512  L1 2.80e-02 1.16e-02 4.41e-03 1.82e-03 | ord  1.27  1.40  1.27
pulse_k12.8  kn03 be       C f1     lev 32-64-128-256-512  L1 7.54e-02 3.17e-02 1.23e-02 5.02e-03 | ord  1.25  1.37  1.30
pulse_k12.8  kn03 hesdirk2 X e      lev 32-64-128-256-512  L1 2.44e-02 6.51e-03 1.35e-03 3.13e-04 | ord  1.90  2.27  2.11
pulse_k12.8  kn03 hesdirk2 X f1     lev 32-64-128-256-512  L1 5.69e-02 1.68e-02 3.64e-03 6.98e-04 | ord  1.76  2.20  2.38
pulse_k12.8  kn03 hesdirk2 T e      lev 1-2-4-8-16  L1 3.97e-03 1.93e-03 9.53e-04 4.73e-04 | ord  1.04  1.02  1.01
pulse_k12.8  kn03 hesdirk2 T f1     lev 1-2-4-8-16  L1 8.73e-03 4.30e-03 2.13e-03 1.06e-03 | ord  1.02  1.01  1.01
pulse_k12.8  kn03 hesdirk2 C e      lev 32-64-128-256-512  L1 2.36e-02 7.36e-03 1.88e-03 5.41e-04 | ord  1.68  1.97  1.80
pulse_k12.8  kn03 hesdirk2 C f1     lev 32-64-128-256-512  L1 6.12e-02 1.98e-02 5.28e-03 1.55e-03 | ord  1.63  1.90  1.77
pulse_k12.8  kn1  be       X e      lev 32-64-128-256-512  L1 2.36e-02 6.38e-03 1.35e-03 3.40e-04 | ord  1.89  2.24  1.99
pulse_k12.8  kn1  be       X f1     lev 32-64-128-256-512  L1 5.56e-02 1.63e-02 3.61e-03 8.72e-04 | ord  1.77  2.17  2.05
pulse_k12.8  kn1  be       T e      lev 1-2-4-8-16  L1 1.36e-02 7.26e-03 3.76e-03 1.92e-03 | ord  0.91  0.95  0.97
pulse_k12.8  kn1  be       T f1     lev 1-2-4-8-16  L1 3.20e-02 1.76e-02 9.27e-03 4.76e-03 | ord  0.86  0.93  0.96
pulse_k12.8  kn1  be       C e      lev 32-64-128-256-512  L1 2.81e-02 1.17e-02 4.51e-03 2.04e-03 | ord  1.27  1.37  1.15
pulse_k12.8  kn1  be       C f1     lev 32-64-128-256-512  L1 7.55e-02 3.19e-02 1.25e-02 5.48e-03 | ord  1.24  1.34  1.19
pulse_k12.8  kn1  hesdirk2 X e      lev 32-64-128-256-512  L1 2.44e-02 6.51e-03 1.38e-03 3.09e-04 | ord  1.90  2.24  2.16
pulse_k12.8  kn1  hesdirk2 X f1     lev 32-64-128-256-512  L1 5.69e-02 1.68e-02 3.73e-03 8.25e-04 | ord  1.76  2.17  2.18
pulse_k12.8  kn1  hesdirk2 T e      lev 1-2-4-8-16  L1 3.80e-03 1.84e-03 9.01e-04 4.46e-04 | ord  1.05  1.03  1.01
pulse_k12.8  kn1  hesdirk2 T f1     lev 1-2-4-8-16  L1 8.27e-03 4.05e-03 2.00e-03 9.94e-04 | ord  1.03  1.02  1.01
pulse_k12.8  kn1  hesdirk2 C e      lev 32-64-128-256-512  L1 2.36e-02 7.37e-03 1.91e-03 6.03e-04 | ord  1.68  1.95  1.67
pulse_k12.8  kn1  hesdirk2 C f1     lev 32-64-128-256-512  L1 6.12e-02 1.98e-02 5.35e-03 1.68e-03 | ord  1.63  1.89  1.67
pulse_k12.8  knp  be       X e      lev 32-64-128-256-512  L1 2.36e-02 6.34e-03 1.26e-03 6.84e-04 | ord  1.89  2.33  0.88
pulse_k12.8  knp  be       X f1     lev 32-64-128-256-512  L1 5.56e-02 1.62e-02 3.28e-03 1.85e-03 | ord  1.78  2.30  0.83
pulse_k12.8  knp  be       T e      lev 1-2-4-8-16  L1 1.42e-02 7.65e-03 3.97e-03 2.03e-03 | ord  0.90  0.94  0.97
pulse_k12.8  knp  be       T f1     lev 1-2-4-8-16  L1 3.33e-02 1.84e-02 9.73e-03 5.00e-03 | ord  0.85  0.92  0.96
pulse_k12.8  knp  be       C e      lev 32-64-128-256-512  L1 2.80e-02 1.16e-02 4.41e-03 1.82e-03 | ord  1.27  1.40  1.27
pulse_k12.8  knp  be       C f1     lev 32-64-128-256-512  L1 7.54e-02 3.17e-02 1.23e-02 5.01e-03 | ord  1.25  1.37  1.30
pulse_k12.8  knp  hesdirk2 X e      lev 32-64-128-256-512  L1 2.44e-02 6.50e-03 1.35e-03 3.09e-04 | ord  1.90  2.27  2.12
pulse_k12.8  knp  hesdirk2 X f1     lev 32-64-128-256-512  L1 5.69e-02 1.67e-02 3.63e-03 6.81e-04 | ord  1.76  2.21  2.41
pulse_k12.8  knp  hesdirk2 T e      lev 1-2-4-8-16  L1 4.01e-03 1.95e-03 9.63e-04 4.78e-04 | ord  1.04  1.02  1.01
pulse_k12.8  knp  hesdirk2 T f1     lev 1-2-4-8-16  L1 8.78e-03 4.33e-03 2.14e-03 1.07e-03 | ord  1.02  1.01  1.01
pulse_k12.8  knp  hesdirk2 C e      lev 32-64-128-256-512  L1 2.36e-02 7.36e-03 1.88e-03 5.41e-04 | ord  1.68  1.97  1.80
pulse_k12.8  knp  hesdirk2 C f1     lev 32-64-128-256-512  L1 6.12e-02 1.98e-02 5.28e-03 1.55e-03 | ord  1.63  1.90  1.77
pulse_k1280  cen  be       X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  cen  be       X f1     lev 32-64-128-256-512  L1 2.32e-02 6.07e-03 1.57e-03 4.15e-04 | ord  1.93  1.95  1.92
pulse_k1280  cen  be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k1280  cen  be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k1280  cen  be       C e      lev 32-64-128-256-512  L1 1.57e-02 4.61e-03 1.35e-03 4.36e-04 | ord  1.77  1.78  1.63
pulse_k1280  cen  be       C f1     lev 32-64-128-256-512  L1 2.21e-02 6.83e-03 2.29e-03 8.82e-04 | ord  1.70  1.58  1.38
pulse_k1280  cen  hesdirk2 X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  cen  hesdirk2 X f1     lev 32-64-128-256-512  L1 2.31e-02 6.07e-03 1.58e-03 4.16e-04 | ord  1.93  1.95  1.92
pulse_k1280  cen  hesdirk2 T e      lev 1-2-4-8-16  L1 4.69e-05 1.21e-05 3.13e-06 8.37e-07 | ord  1.96  1.94  1.90
pulse_k1280  cen  hesdirk2 T f1     lev 1-2-4-8-16  L1 1.67e-04 4.28e-05 1.10e-05 2.86e-06 | ord  1.97  1.96  1.94
pulse_k1280  cen  hesdirk2 C e      lev 32-64-128-256-512  L1 1.46e-02 3.93e-03 1.00e-03 2.63e-04 | ord  1.89  1.97  1.93
pulse_k1280  cen  hesdirk2 C f1     lev 32-64-128-256-512  L1 2.30e-02 6.07e-03 1.58e-03 4.17e-04 | ord  1.92  1.94  1.92
pulse_k1280  hr   be       X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 1.01e-03 2.92e-04 | ord  1.89  1.95  1.80
pulse_k1280  hr   be       X f1     lev 32-64-128-256-512  L1 2.32e-02 6.08e-03 1.60e-03 4.69e-04 | ord  1.93  1.93  1.77
pulse_k1280  hr   be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.50e-04 1.75e-04 | ord  0.99  1.00  1.00
pulse_k1280  hr   be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.30e-03 6.52e-04 | ord  0.99  0.99  1.00
pulse_k1280  hr   be       C e      lev 32-64-128-256-512  L1 1.57e-02 4.61e-03 1.36e-03 4.60e-04 | ord  1.76  1.76  1.56
pulse_k1280  hr   be       C f1     lev 32-64-128-256-512  L1 2.22e-02 6.84e-03 2.31e-03 9.23e-04 | ord  1.70  1.57  1.32
pulse_k1280  hr   hesdirk2 X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 1.01e-03 2.89e-04 | ord  1.89  1.95  1.81
pulse_k1280  hr   hesdirk2 X f1     lev 32-64-128-256-512  L1 2.31e-02 6.08e-03 1.60e-03 4.64e-04 | ord  1.93  1.93  1.78
pulse_k1280  hr   hesdirk2 T e      lev 1-2-4-8-16  L1 4.62e-05 1.10e-05 1.87e-06 3.67e-06 | ord  2.08  2.55 -0.97
pulse_k1280  hr   hesdirk2 T f1     lev 1-2-4-8-16  L1 1.66e-04 4.14e-05 1.23e-05 1.36e-05 | ord  2.01  1.75 -0.14
pulse_k1280  hr   hesdirk2 C e      lev 32-64-128-256-512  L1 1.46e-02 3.94e-03 1.01e-03 2.82e-04 | ord  1.89  1.96  1.84
pulse_k1280  hr   hesdirk2 C f1     lev 32-64-128-256-512  L1 2.30e-02 6.09e-03 1.60e-03 4.50e-04 | ord  1.92  1.93  1.83
pulse_k1280  sq   be       X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  sq   be       X f1     lev 32-64-128-256-512  L1 2.32e-02 6.07e-03 1.57e-03 4.15e-04 | ord  1.93  1.95  1.92
pulse_k1280  sq   be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k1280  sq   be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k1280  sq   be       C e      lev 32-64-128-256-512  L1 1.57e-02 4.61e-03 1.35e-03 4.36e-04 | ord  1.77  1.78  1.63
pulse_k1280  sq   be       C f1     lev 32-64-128-256-512  L1 2.21e-02 6.83e-03 2.29e-03 8.82e-04 | ord  1.70  1.58  1.38
pulse_k1280  sq   hesdirk2 X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  sq   hesdirk2 X f1     lev 32-64-128-256-512  L1 2.31e-02 6.07e-03 1.58e-03 4.16e-04 | ord  1.93  1.95  1.92
pulse_k1280  sq   hesdirk2 T e      lev 1-2-4-8-16  L1 4.69e-05 1.21e-05 3.13e-06 8.37e-07 | ord  1.96  1.94  1.90
pulse_k1280  sq   hesdirk2 T f1     lev 1-2-4-8-16  L1 1.67e-04 4.28e-05 1.10e-05 2.86e-06 | ord  1.97  1.96  1.94
pulse_k1280  sq   hesdirk2 C e      lev 32-64-128-256-512  L1 1.46e-02 3.93e-03 1.00e-03 2.63e-04 | ord  1.89  1.97  1.93
pulse_k1280  sq   hesdirk2 C f1     lev 32-64-128-256-512  L1 2.30e-02 6.07e-03 1.58e-03 4.17e-04 | ord  1.92  1.94  1.92
pulse_k1280  sqf  be       X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  sqf  be       X f1     lev 32-64-128-256-512  L1 2.32e-02 6.07e-03 1.57e-03 4.15e-04 | ord  1.93  1.95  1.92
pulse_k1280  sqf  be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k1280  sqf  be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k1280  sqf  be       C e      lev 32-64-128-256-512  L1 1.57e-02 4.61e-03 1.35e-03 4.36e-04 | ord  1.77  1.78  1.63
pulse_k1280  sqf  be       C f1     lev 32-64-128-256-512  L1 2.21e-02 6.83e-03 2.29e-03 8.82e-04 | ord  1.70  1.58  1.38
pulse_k1280  sqf  hesdirk2 X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  sqf  hesdirk2 X f1     lev 32-64-128-256-512  L1 2.31e-02 6.07e-03 1.58e-03 4.16e-04 | ord  1.93  1.95  1.92
pulse_k1280  sqf  hesdirk2 T e      lev 1-2-4-8-16  L1 4.69e-05 1.21e-05 3.13e-06 8.37e-07 | ord  1.96  1.94  1.90
pulse_k1280  sqf  hesdirk2 T f1     lev 1-2-4-8-16  L1 1.67e-04 4.28e-05 1.10e-05 2.86e-06 | ord  1.97  1.96  1.94
pulse_k1280  sqf  hesdirk2 C e      lev 32-64-128-256-512  L1 1.46e-02 3.93e-03 1.00e-03 2.63e-04 | ord  1.89  1.97  1.93
pulse_k1280  sqf  hesdirk2 C f1     lev 32-64-128-256-512  L1 2.30e-02 6.07e-03 1.58e-03 4.17e-04 | ord  1.92  1.94  1.92
pulse_k1280  sqv  be       X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  sqv  be       X f1     lev 32-64-128-256-512  L1 2.32e-02 6.07e-03 1.57e-03 4.15e-04 | ord  1.93  1.95  1.92
pulse_k1280  sqv  be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k1280  sqv  be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k1280  sqv  be       C e      lev 32-64-128-256-512  L1 1.57e-02 4.61e-03 1.35e-03 4.36e-04 | ord  1.77  1.78  1.63
pulse_k1280  sqv  be       C f1     lev 32-64-128-256-512  L1 2.21e-02 6.83e-03 2.29e-03 8.82e-04 | ord  1.70  1.58  1.38
pulse_k1280  sqv  hesdirk2 X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  sqv  hesdirk2 X f1     lev 32-64-128-256-512  L1 2.31e-02 6.07e-03 1.58e-03 4.16e-04 | ord  1.93  1.95  1.92
pulse_k1280  sqv  hesdirk2 T e      lev 1-2-4-8-16  L1 4.69e-05 1.21e-05 3.13e-06 8.37e-07 | ord  1.96  1.94  1.90
pulse_k1280  sqv  hesdirk2 T f1     lev 1-2-4-8-16  L1 1.67e-04 4.28e-05 1.10e-05 2.86e-06 | ord  1.97  1.96  1.94
pulse_k1280  sqv  hesdirk2 C e      lev 32-64-128-256-512  L1 1.46e-02 3.93e-03 1.00e-03 2.63e-04 | ord  1.89  1.97  1.93
pulse_k1280  sqv  hesdirk2 C f1     lev 32-64-128-256-512  L1 2.30e-02 6.07e-03 1.58e-03 4.17e-04 | ord  1.92  1.94  1.92
pulse_k1280  sqa  be       X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  sqa  be       X f1     lev 32-64-128-256-512  L1 2.32e-02 6.07e-03 1.57e-03 4.15e-04 | ord  1.93  1.95  1.92
pulse_k1280  sqa  be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k1280  sqa  be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k1280  sqa  be       C e      lev 32-64-128-256-512  L1 1.57e-02 4.61e-03 1.35e-03 4.36e-04 | ord  1.77  1.78  1.63
pulse_k1280  sqa  be       C f1     lev 32-64-128-256-512  L1 2.21e-02 6.83e-03 2.29e-03 8.82e-04 | ord  1.70  1.58  1.38
pulse_k1280  sqa  hesdirk2 X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  sqa  hesdirk2 X f1     lev 32-64-128-256-512  L1 2.31e-02 6.07e-03 1.58e-03 4.16e-04 | ord  1.93  1.95  1.92
pulse_k1280  sqa  hesdirk2 T e      lev 1-2-4-8-16  L1 4.69e-05 1.21e-05 3.13e-06 8.37e-07 | ord  1.96  1.94  1.90
pulse_k1280  sqa  hesdirk2 T f1     lev 1-2-4-8-16  L1 1.67e-04 4.28e-05 1.10e-05 2.86e-06 | ord  1.97  1.96  1.94
pulse_k1280  sqa  hesdirk2 C e      lev 32-64-128-256-512  L1 1.46e-02 3.93e-03 1.00e-03 2.63e-04 | ord  1.89  1.97  1.93
pulse_k1280  sqa  hesdirk2 C f1     lev 32-64-128-256-512  L1 2.30e-02 6.07e-03 1.58e-03 4.17e-04 | ord  1.92  1.94  1.92
pulse_k1280  aq   be       X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  aq   be       X f1     lev 32-64-128-256-512  L1 2.32e-02 6.07e-03 1.57e-03 4.15e-04 | ord  1.93  1.95  1.92
pulse_k1280  aq   be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k1280  aq   be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k1280  aq   be       C e      lev 32-64-128-256-512  L1 1.57e-02 4.61e-03 1.35e-03 4.36e-04 | ord  1.77  1.78  1.63
pulse_k1280  aq   be       C f1     lev 32-64-128-256-512  L1 2.21e-02 6.83e-03 2.29e-03 8.82e-04 | ord  1.70  1.58  1.38
pulse_k1280  aq   hesdirk2 X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  aq   hesdirk2 X f1     lev 32-64-128-256-512  L1 2.31e-02 6.07e-03 1.58e-03 4.16e-04 | ord  1.93  1.95  1.92
pulse_k1280  aq   hesdirk2 T e      lev 1-2-4-8-16  L1 4.68e-05 1.20e-05 3.10e-06 8.13e-07 | ord  1.96  1.95  1.93
pulse_k1280  aq   hesdirk2 T f1     lev 1-2-4-8-16  L1 1.67e-04 4.27e-05 1.09e-05 2.82e-06 | ord  1.97  1.97  1.95
pulse_k1280  aq   hesdirk2 C e      lev 32-64-128-256-512  L1 1.46e-02 3.93e-03 1.00e-03 2.63e-04 | ord  1.89  1.97  1.93
pulse_k1280  aq   hesdirk2 C f1     lev 32-64-128-256-512  L1 2.30e-02 6.07e-03 1.58e-03 4.16e-04 | ord  1.92  1.94  1.92
pulse_k1280  kn01 be       X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  kn01 be       X f1     lev 32-64-128-256-512  L1 2.32e-02 6.07e-03 1.57e-03 4.15e-04 | ord  1.93  1.95  1.92
pulse_k1280  kn01 be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k1280  kn01 be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k1280  kn01 be       C e      lev 32-64-128-256-512  L1 1.57e-02 4.61e-03 1.35e-03 4.36e-04 | ord  1.77  1.78  1.63
pulse_k1280  kn01 be       C f1     lev 32-64-128-256-512  L1 2.21e-02 6.83e-03 2.29e-03 8.82e-04 | ord  1.70  1.58  1.38
pulse_k1280  kn01 hesdirk2 X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  kn01 hesdirk2 X f1     lev 32-64-128-256-512  L1 2.31e-02 6.07e-03 1.58e-03 4.16e-04 | ord  1.93  1.95  1.92
pulse_k1280  kn01 hesdirk2 T e      lev 1-2-4-8-16  L1 4.69e-05 1.21e-05 3.13e-06 8.37e-07 | ord  1.96  1.94  1.90
pulse_k1280  kn01 hesdirk2 T f1     lev 1-2-4-8-16  L1 1.67e-04 4.28e-05 1.10e-05 2.86e-06 | ord  1.97  1.96  1.94
pulse_k1280  kn01 hesdirk2 C e      lev 32-64-128-256-512  L1 1.46e-02 3.93e-03 1.00e-03 2.63e-04 | ord  1.89  1.97  1.93
pulse_k1280  kn01 hesdirk2 C f1     lev 32-64-128-256-512  L1 2.30e-02 6.07e-03 1.58e-03 4.17e-04 | ord  1.92  1.94  1.92
pulse_k1280  kn03 be       X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  kn03 be       X f1     lev 32-64-128-256-512  L1 2.32e-02 6.07e-03 1.57e-03 4.15e-04 | ord  1.93  1.95  1.92
pulse_k1280  kn03 be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k1280  kn03 be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k1280  kn03 be       C e      lev 32-64-128-256-512  L1 1.57e-02 4.61e-03 1.35e-03 4.36e-04 | ord  1.77  1.78  1.63
pulse_k1280  kn03 be       C f1     lev 32-64-128-256-512  L1 2.21e-02 6.83e-03 2.29e-03 8.82e-04 | ord  1.70  1.58  1.38
pulse_k1280  kn03 hesdirk2 X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  kn03 hesdirk2 X f1     lev 32-64-128-256-512  L1 2.31e-02 6.07e-03 1.58e-03 4.16e-04 | ord  1.93  1.95  1.92
pulse_k1280  kn03 hesdirk2 T e      lev 1-2-4-8-16  L1 4.69e-05 1.21e-05 3.13e-06 8.37e-07 | ord  1.96  1.94  1.90
pulse_k1280  kn03 hesdirk2 T f1     lev 1-2-4-8-16  L1 1.67e-04 4.28e-05 1.10e-05 2.86e-06 | ord  1.97  1.96  1.94
pulse_k1280  kn03 hesdirk2 C e      lev 32-64-128-256-512  L1 1.46e-02 3.93e-03 1.00e-03 2.63e-04 | ord  1.89  1.97  1.93
pulse_k1280  kn03 hesdirk2 C f1     lev 32-64-128-256-512  L1 2.30e-02 6.07e-03 1.58e-03 4.17e-04 | ord  1.92  1.94  1.92
pulse_k1280  kn1  be       X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  kn1  be       X f1     lev 32-64-128-256-512  L1 2.32e-02 6.07e-03 1.57e-03 4.15e-04 | ord  1.93  1.95  1.92
pulse_k1280  kn1  be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k1280  kn1  be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k1280  kn1  be       C e      lev 32-64-128-256-512  L1 1.57e-02 4.61e-03 1.35e-03 4.36e-04 | ord  1.77  1.78  1.63
pulse_k1280  kn1  be       C f1     lev 32-64-128-256-512  L1 2.21e-02 6.83e-03 2.29e-03 8.82e-04 | ord  1.70  1.58  1.38
pulse_k1280  kn1  hesdirk2 X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  kn1  hesdirk2 X f1     lev 32-64-128-256-512  L1 2.31e-02 6.07e-03 1.58e-03 4.16e-04 | ord  1.93  1.95  1.92
pulse_k1280  kn1  hesdirk2 T e      lev 1-2-4-8-16  L1 4.69e-05 1.21e-05 3.13e-06 8.37e-07 | ord  1.96  1.94  1.90
pulse_k1280  kn1  hesdirk2 T f1     lev 1-2-4-8-16  L1 1.67e-04 4.28e-05 1.10e-05 2.86e-06 | ord  1.97  1.96  1.94
pulse_k1280  kn1  hesdirk2 C e      lev 32-64-128-256-512  L1 1.46e-02 3.93e-03 1.00e-03 2.63e-04 | ord  1.89  1.97  1.93
pulse_k1280  kn1  hesdirk2 C f1     lev 32-64-128-256-512  L1 2.30e-02 6.07e-03 1.58e-03 4.17e-04 | ord  1.92  1.94  1.92
pulse_k1280  knp  be       X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  knp  be       X f1     lev 32-64-128-256-512  L1 2.32e-02 6.07e-03 1.57e-03 4.15e-04 | ord  1.93  1.95  1.92
pulse_k1280  knp  be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k1280  knp  be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k1280  knp  be       C e      lev 32-64-128-256-512  L1 1.57e-02 4.61e-03 1.35e-03 4.36e-04 | ord  1.77  1.78  1.63
pulse_k1280  knp  be       C f1     lev 32-64-128-256-512  L1 2.21e-02 6.83e-03 2.29e-03 8.82e-04 | ord  1.70  1.58  1.38
pulse_k1280  knp  hesdirk2 X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  knp  hesdirk2 X f1     lev 32-64-128-256-512  L1 2.31e-02 6.07e-03 1.58e-03 4.16e-04 | ord  1.93  1.95  1.92
pulse_k1280  knp  hesdirk2 T e      lev 1-2-4-8-16  L1 4.69e-05 1.21e-05 3.13e-06 8.37e-07 | ord  1.96  1.94  1.90
pulse_k1280  knp  hesdirk2 T f1     lev 1-2-4-8-16  L1 1.67e-04 4.28e-05 1.10e-05 2.86e-06 | ord  1.97  1.96  1.94
pulse_k1280  knp  hesdirk2 C e      lev 32-64-128-256-512  L1 1.46e-02 3.93e-03 1.00e-03 2.63e-04 | ord  1.89  1.97  1.93
pulse_k1280  knp  hesdirk2 C f1     lev 32-64-128-256-512  L1 2.30e-02 6.07e-03 1.58e-03 4.17e-04 | ord  1.92  1.94  1.92
rw_t10       cen  be       X e      lev 32-64-128-256-512  L1 2.28e-02 5.07e-03 1.23e-03 3.17e-04 | ord  2.17  2.04  1.96
rw_t10       cen  be       X f1     lev 32-64-128-256-512  L1 2.01e-01 1.48e-01 1.37e-01 1.10e-03 | ord  0.44  0.12  6.95
rw_t10       cen  be       X dens   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 2.98e-04 | ord  2.18  2.14  2.16
rw_t10       cen  be       X velx   lev 32-64-128-256-512  L1 2.72e-02 6.07e-03 1.37e-03 3.13e-04 | ord  2.17  2.15  2.13
rw_t10       cen  be       X eint   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 2.98e-04 | ord  2.18  2.14  2.16
rw_t10       cen  be       T e      lev 1-2-4-8-16  L1 7.83e-02 4.28e-02 2.24e-02 1.15e-02 | ord  0.87  0.93  0.96
rw_t10       cen  be       T f1     lev 1-2-4-8-16  L1 6.62e-02 3.55e-02 1.85e-02 9.38e-03 | ord  0.90  0.94  0.98
rw_t10       cen  be       T dens   lev 1-2-4-8-16  L1 6.75e-03 3.37e-03 1.68e-03 8.41e-04 | ord  1.00  1.00  1.00
rw_t10       cen  be       T velx   lev 1-2-4-8-16  L1 6.41e-03 3.20e-03 1.60e-03 8.00e-04 | ord  1.00  1.00  1.00
rw_t10       cen  be       T eint   lev 1-2-4-8-16  L1 6.78e-03 3.39e-03 1.69e-03 8.45e-04 | ord  1.00  1.00  1.00
rw_t10       cen  be       C e      lev 32-64-128-256-512  L1 2.00e-01 1.33e-01 7.95e-02 4.36e-02 | ord  0.59  0.74  0.87
rw_t10       cen  be       C f1     lev 32-64-128-256-512  L1 2.11e-01 1.22e-01 6.67e-02 3.63e-02 | ord  0.79  0.87  0.88
rw_t10       cen  be       C dens   lev 32-64-128-256-512  L1 4.40e-02 1.51e-02 6.85e-03 3.40e-03 | ord  1.54  1.14  1.01
rw_t10       cen  be       C velx   lev 32-64-128-256-512  L1 4.01e-02 1.37e-02 6.48e-03 3.22e-03 | ord  1.55  1.08  1.01
rw_t10       cen  be       C eint   lev 32-64-128-256-512  L1 4.45e-02 1.54e-02 6.91e-03 3.42e-03 | ord  1.54  1.15  1.01
rw_t10       cen  hesdirk2 X e      lev 32-64-128-256-512  L1 1.54e-02 4.10e-03 1.05e-03 2.81e-04 | ord  1.91  1.97  1.90
rw_t10       cen  hesdirk2 X f1     lev 32-64-128-256-512  L1 1.67e-01 4.76e-02 4.50e-02 1.10e-03 | ord  1.81  0.08  5.36
rw_t10       cen  hesdirk2 X dens   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.97e-04 | ord  2.18  2.13  2.19
rw_t10       cen  hesdirk2 X velx   lev 32-64-128-256-512  L1 2.74e-02 6.13e-03 1.39e-03 3.16e-04 | ord  2.16  2.14  2.14
rw_t10       cen  hesdirk2 X eint   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.97e-04 | ord  2.18  2.13  2.19
rw_t10       cen  hesdirk2 T e      lev 1-2-4-8-16  L1 8.57e-04 4.25e-04 2.12e-04 1.06e-04 | ord  1.01  1.00  1.00
rw_t10       cen  hesdirk2 T f1     lev 1-2-4-8-16  L1 7.30e-04 2.67e-04 6.35e-04 8.44e-04 | ord  1.45 -1.25 -0.41
rw_t10       cen  hesdirk2 T dens   lev 1-2-4-8-16  L1 1.27e-03 6.33e-04 3.16e-04 1.58e-04 | ord  1.00  1.00  1.00
rw_t10       cen  hesdirk2 T velx   lev 1-2-4-8-16  L1 1.29e-03 6.38e-04 3.17e-04 1.58e-04 | ord  1.02  1.01  1.00
rw_t10       cen  hesdirk2 T eint   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.00  1.00  1.00
rw_t10       cen  hesdirk2 C e      lev 32-64-128-256-512  L1 1.43e-02 2.81e-03 5.26e-04 2.42e-04 | ord  2.35  2.42  1.12
rw_t10       cen  hesdirk2 C f1     lev 32-64-128-256-512  L1 1.57e-01 2.81e-02 3.81e-03 1.32e-03 | ord  2.49  2.88  1.52
rw_t10       cen  hesdirk2 C dens   lev 32-64-128-256-512  L1 2.86e-02 6.61e-03 1.70e-03 6.19e-04 | ord  2.12  1.96  1.46
rw_t10       cen  hesdirk2 C velx   lev 32-64-128-256-512  L1 2.87e-02 6.43e-03 1.53e-03 6.17e-04 | ord  2.16  2.07  1.31
rw_t10       cen  hesdirk2 C eint   lev 32-64-128-256-512  L1 2.87e-02 6.61e-03 1.70e-03 6.19e-04 | ord  2.12  1.96  1.46
rw_t10       hr   be       X e      lev 32-64-128-256-512  L1 2.28e-02 3.93e-02 1.14e-01 1.96e-01 | ord -0.79 -1.53 -0.79
rw_t10       hr   be       X f1     lev 32-64-128-256-512  L1 2.06e-01 2.39e-01 5.75e-01 1.02e+00 | ord -0.22 -1.27 -0.83
rw_t10       hr   be       X dens   lev 32-64-128-256-512  L1 2.66e-02 6.24e-03 2.68e-03 2.73e-03 | ord  2.09  1.22 -0.03
rw_t10       hr   be       X velx   lev 32-64-128-256-512  L1 2.73e-02 6.44e-03 3.36e-03 4.84e-03 | ord  2.08  0.94 -0.53
rw_t10       hr   be       X eint   lev 32-64-128-256-512  L1 2.67e-02 6.30e-03 2.86e-03 2.99e-03 | ord  2.08  1.14 -0.06
rw_t10       hr   be       T e      lev 1-2-4-8-16  L1 8.43e-02 3.58e-02 1.34e-02 2.90e-03 | ord  1.23  1.42  2.21
rw_t10       hr   be       T f1     lev 1-2-4-8-16  L1 6.77e-02 5.03e-02 6.64e-02 7.34e-02 | ord  0.43 -0.40 -0.14
rw_t10       hr   be       T dens   lev 1-2-4-8-16  L1 5.03e-03 2.54e-03 1.30e-03 6.84e-04 | ord  0.99  0.96  0.93
rw_t10       hr   be       T velx   lev 1-2-4-8-16  L1 5.22e-03 2.64e-03 1.36e-03 7.36e-04 | ord  0.98  0.96  0.88
rw_t10       hr   be       T eint   lev 1-2-4-8-16  L1 5.05e-03 2.55e-03 1.30e-03 6.84e-04 | ord  0.99  0.97  0.93
rw_t10       hr   be       C e      lev 32-64-128-256-512  L1 2.63e-01 2.63e-01 2.12e-01 3.34e-02 | ord  0.00  0.31  2.67
rw_t10       hr   be       C f1     lev 32-64-128-256-512  L1 2.58e-01 2.85e-01 7.16e-01 1.06e+00 | ord -0.14 -1.33 -0.57
rw_t10       hr   be       C dens   lev 32-64-128-256-512  L1 4.52e-02 1.71e-02 8.08e-03 2.60e-03 | ord  1.40  1.08  1.64
rw_t10       hr   be       C velx   lev 32-64-128-256-512  L1 4.12e-02 1.52e-02 7.71e-03 4.32e-03 | ord  1.44  0.98  0.84
rw_t10       hr   be       C eint   lev 32-64-128-256-512  L1 4.58e-02 1.75e-02 8.23e-03 2.60e-03 | ord  1.39  1.09  1.66
rw_t10       hr   hesdirk2 X e      lev 32-64-128-256-512  L1 1.54e-02 5.30e-03 1.69e-02 3.83e-02 | ord  1.54 -1.67 -1.18
rw_t10       hr   hesdirk2 X f1     lev 32-64-128-256-512  L1 1.68e-01 1.53e-01 4.42e-01 8.66e-01 | ord  0.14 -1.53 -0.97
rw_t10       hr   hesdirk2 X dens   lev 32-64-128-256-512  L1 2.67e-02 5.86e-03 1.22e-03 4.96e-04 | ord  2.19  2.27  1.30
rw_t10       hr   hesdirk2 X velx   lev 32-64-128-256-512  L1 2.74e-02 6.11e-03 2.39e-03 4.60e-03 | ord  2.16  1.36 -0.95
rw_t10       hr   hesdirk2 X eint   lev 32-64-128-256-512  L1 2.67e-02 5.86e-03 1.20e-03 5.21e-04 | ord  2.19  2.29  1.20
rw_t10       hr   hesdirk2 T e      lev 1-2-4-8-16  L1 6.87e-04 1.25e-03 1.31e-03 1.04e-03 | ord -0.86 -0.07  0.34
rw_t10       hr   hesdirk2 T f1     lev 1-2-4-8-16  L1 6.15e-02 7.27e-02 6.83e-02 5.18e-02 | ord -0.24  0.09  0.40
rw_t10       hr   hesdirk2 T dens   lev 1-2-4-8-16  L1 1.28e-03 6.48e-04 3.31e-04 1.69e-04 | ord  0.99  0.97  0.97
rw_t10       hr   hesdirk2 T velx   lev 1-2-4-8-16  L1 1.31e-03 6.93e-04 4.10e-04 2.71e-04 | ord  0.92  0.76  0.60
rw_t10       hr   hesdirk2 T eint   lev 1-2-4-8-16  L1 1.28e-03 6.48e-04 3.31e-04 1.69e-04 | ord  0.99  0.97  0.97
rw_t10       hr   hesdirk2 C e      lev 32-64-128-256-512  L1 1.46e-02 8.81e-03 2.06e-02 3.22e-02 | ord  0.73 -1.22 -0.65
rw_t10       hr   hesdirk2 C f1     lev 32-64-128-256-512  L1 1.79e-01 2.28e-01 5.01e-01 6.77e-01 | ord -0.35 -1.14 -0.43
rw_t10       hr   hesdirk2 C dens   lev 32-64-128-256-512  L1 2.86e-02 6.53e-03 1.43e-03 2.77e-04 | ord  2.13  2.19  2.36
rw_t10       hr   hesdirk2 C velx   lev 32-64-128-256-512  L1 2.87e-02 6.40e-03 2.78e-03 3.61e-03 | ord  2.16  1.21 -0.38
rw_t10       hr   hesdirk2 C eint   lev 32-64-128-256-512  L1 2.87e-02 6.52e-03 1.40e-03 3.08e-04 | ord  2.13  2.22  2.18
rw_t10       sq   be       X e      lev 32-64-128-256-512  L1 2.28e-02 5.07e-03 1.25e-03 1.14e-02 | ord  2.17  2.02 -3.19
rw_t10       sq   be       X f1     lev 32-64-128-256-512  L1 2.01e-01 1.48e-01 1.37e-01 7.04e-02 | ord  0.44  0.11  0.96
rw_t10       sq   be       X dens   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 4.32e-04 | ord  2.18  2.14  1.62
rw_t10       sq   be       X velx   lev 32-64-128-256-512  L1 2.72e-02 6.07e-03 1.37e-03 4.65e-04 | ord  2.17  2.15  1.56
rw_t10       sq   be       X eint   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 4.52e-04 | ord  2.18  2.14  1.56
rw_t10       sq   be       T e      lev 1-2-4-8-16  L1 7.83e-02 4.28e-02 2.24e-02 1.15e-02 | ord  0.87  0.93  0.96
rw_t10       sq   be       T f1     lev 1-2-4-8-16  L1 6.68e-02 3.57e-02 1.85e-02 9.37e-03 | ord  0.90  0.95  0.98
rw_t10       sq   be       T dens   lev 1-2-4-8-16  L1 6.76e-03 3.37e-03 1.68e-03 8.41e-04 | ord  1.00  1.00  1.00
rw_t10       sq   be       T velx   lev 1-2-4-8-16  L1 6.41e-03 3.20e-03 1.60e-03 8.01e-04 | ord  1.00  1.00  1.00
rw_t10       sq   be       T eint   lev 1-2-4-8-16  L1 6.79e-03 3.39e-03 1.69e-03 8.45e-04 | ord  1.00  1.00  1.00
rw_t10       sq   be       C e      lev 32-64-128-256-512  L1 2.03e-01 1.38e-01 7.98e-02 5.34e-02 | ord  0.56  0.79  0.58
rw_t10       sq   be       C f1     lev 32-64-128-256-512  L1 2.09e-01 1.24e-01 6.77e-02 6.49e-02 | ord  0.76  0.87  0.06
rw_t10       sq   be       C dens   lev 32-64-128-256-512  L1 4.39e-02 1.52e-02 6.87e-03 3.53e-03 | ord  1.53  1.15  0.96
rw_t10       sq   be       C velx   lev 32-64-128-256-512  L1 4.01e-02 1.37e-02 6.48e-03 3.33e-03 | ord  1.55  1.08  0.96
rw_t10       sq   be       C eint   lev 32-64-128-256-512  L1 4.44e-02 1.55e-02 6.93e-03 3.56e-03 | ord  1.52  1.16  0.96
rw_t10       sq   hesdirk2 X e      lev 32-64-128-256-512  L1 1.54e-02 4.10e-03 1.05e-03 2.85e-04 | ord  1.91  1.97  1.88
rw_t10       sq   hesdirk2 X f1     lev 32-64-128-256-512  L1 1.68e-01 4.76e-02 4.51e-02 2.29e-03 | ord  1.82  0.08  4.30
rw_t10       sq   hesdirk2 X dens   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.93e-04 | ord  2.18  2.13  2.21
rw_t10       sq   hesdirk2 X velx   lev 32-64-128-256-512  L1 2.74e-02 6.13e-03 1.39e-03 3.12e-04 | ord  2.16  2.14  2.16
rw_t10       sq   hesdirk2 X eint   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.93e-04 | ord  2.18  2.13  2.21
rw_t10       sq   hesdirk2 T e      lev 1-2-4-8-16  L1 8.55e-04 4.27e-04 4.75e-04 5.20e-04 | ord  1.00 -0.15 -0.13
rw_t10       sq   hesdirk2 T f1     lev 1-2-4-8-16  L1 4.40e-04 4.94e-04 2.33e-04 8.51e-05 | ord -0.17  1.08  1.45
rw_t10       sq   hesdirk2 T dens   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       sq   hesdirk2 T velx   lev 1-2-4-8-16  L1 1.29e-03 6.37e-04 3.17e-04 1.58e-04 | ord  1.02  1.01  1.01
rw_t10       sq   hesdirk2 T eint   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       sq   hesdirk2 C e      lev 32-64-128-256-512  L1 1.43e-02 2.81e-03 5.28e-04 3.16e-04 | ord  2.35  2.41  0.74
rw_t10       sq   hesdirk2 C f1     lev 32-64-128-256-512  L1 1.57e-01 2.82e-02 4.02e-03 2.24e-03 | ord  2.48  2.81  0.84
rw_t10       sq   hesdirk2 C dens   lev 32-64-128-256-512  L1 2.86e-02 6.61e-03 1.70e-03 6.13e-04 | ord  2.12  1.96  1.47
rw_t10       sq   hesdirk2 C velx   lev 32-64-128-256-512  L1 2.87e-02 6.43e-03 1.53e-03 6.16e-04 | ord  2.16  2.07  1.31
rw_t10       sq   hesdirk2 C eint   lev 32-64-128-256-512  L1 2.87e-02 6.61e-03 1.70e-03 6.13e-04 | ord  2.12  1.96  1.47
rw_t10       sqf  be       X e      lev 32-64-128-256-512  L1 2.28e-02 5.07e-03 1.25e-03 1.14e-02 | ord  2.17  2.02 -3.19
rw_t10       sqf  be       X f1     lev 32-64-128-256-512  L1 2.01e-01 1.48e-01 1.37e-01 7.04e-02 | ord  0.44  0.11  0.96
rw_t10       sqf  be       X dens   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 4.32e-04 | ord  2.18  2.14  1.62
rw_t10       sqf  be       X velx   lev 32-64-128-256-512  L1 2.72e-02 6.07e-03 1.37e-03 4.65e-04 | ord  2.17  2.15  1.56
rw_t10       sqf  be       X eint   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 4.52e-04 | ord  2.18  2.14  1.56
rw_t10       sqf  be       T e      lev 1-2-4-8-16  L1 7.83e-02 4.28e-02 2.24e-02 1.15e-02 | ord  0.87  0.93  0.96
rw_t10       sqf  be       T f1     lev 1-2-4-8-16  L1 6.68e-02 3.57e-02 1.85e-02 9.37e-03 | ord  0.90  0.95  0.98
rw_t10       sqf  be       T dens   lev 1-2-4-8-16  L1 6.76e-03 3.37e-03 1.68e-03 8.41e-04 | ord  1.00  1.00  1.00
rw_t10       sqf  be       T velx   lev 1-2-4-8-16  L1 6.41e-03 3.20e-03 1.60e-03 8.01e-04 | ord  1.00  1.00  1.00
rw_t10       sqf  be       T eint   lev 1-2-4-8-16  L1 6.79e-03 3.39e-03 1.69e-03 8.45e-04 | ord  1.00  1.00  1.00
rw_t10       sqf  be       C e      lev 32-64-128-256-512  L1 2.03e-01 1.38e-01 7.98e-02 5.34e-02 | ord  0.56  0.79  0.58
rw_t10       sqf  be       C f1     lev 32-64-128-256-512  L1 2.09e-01 1.24e-01 6.77e-02 6.49e-02 | ord  0.76  0.87  0.06
rw_t10       sqf  be       C dens   lev 32-64-128-256-512  L1 4.39e-02 1.52e-02 6.87e-03 3.53e-03 | ord  1.53  1.15  0.96
rw_t10       sqf  be       C velx   lev 32-64-128-256-512  L1 4.01e-02 1.37e-02 6.48e-03 3.33e-03 | ord  1.55  1.08  0.96
rw_t10       sqf  be       C eint   lev 32-64-128-256-512  L1 4.44e-02 1.55e-02 6.93e-03 3.56e-03 | ord  1.52  1.16  0.96
rw_t10       sqf  hesdirk2 X e      lev 32-64-128-256-512  L1 1.54e-02 4.10e-03 1.05e-03 2.85e-04 | ord  1.91  1.97  1.88
rw_t10       sqf  hesdirk2 X f1     lev 32-64-128-256-512  L1 1.68e-01 4.76e-02 4.51e-02 2.29e-03 | ord  1.82  0.08  4.30
rw_t10       sqf  hesdirk2 X dens   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.93e-04 | ord  2.18  2.13  2.21
rw_t10       sqf  hesdirk2 X velx   lev 32-64-128-256-512  L1 2.74e-02 6.13e-03 1.39e-03 3.12e-04 | ord  2.16  2.14  2.16
rw_t10       sqf  hesdirk2 X eint   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.93e-04 | ord  2.18  2.13  2.21
rw_t10       sqf  hesdirk2 T e      lev 1-2-4-8-16  L1 8.55e-04 4.27e-04 4.75e-04 5.20e-04 | ord  1.00 -0.15 -0.13
rw_t10       sqf  hesdirk2 T f1     lev 1-2-4-8-16  L1 4.40e-04 4.94e-04 2.33e-04 8.51e-05 | ord -0.17  1.08  1.45
rw_t10       sqf  hesdirk2 T dens   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       sqf  hesdirk2 T velx   lev 1-2-4-8-16  L1 1.29e-03 6.37e-04 3.17e-04 1.58e-04 | ord  1.02  1.01  1.01
rw_t10       sqf  hesdirk2 T eint   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       sqf  hesdirk2 C e      lev 32-64-128-256-512  L1 1.43e-02 2.81e-03 5.28e-04 3.16e-04 | ord  2.35  2.41  0.74
rw_t10       sqf  hesdirk2 C f1     lev 32-64-128-256-512  L1 1.57e-01 2.82e-02 4.02e-03 2.24e-03 | ord  2.48  2.81  0.84
rw_t10       sqf  hesdirk2 C dens   lev 32-64-128-256-512  L1 2.86e-02 6.61e-03 1.70e-03 6.13e-04 | ord  2.12  1.96  1.47
rw_t10       sqf  hesdirk2 C velx   lev 32-64-128-256-512  L1 2.87e-02 6.43e-03 1.53e-03 6.16e-04 | ord  2.16  2.07  1.31
rw_t10       sqf  hesdirk2 C eint   lev 32-64-128-256-512  L1 2.87e-02 6.61e-03 1.70e-03 6.13e-04 | ord  2.12  1.96  1.47
rw_t10       sqv  be       X e      lev 32-64-128-256-512  L1 2.28e-02 5.07e-03 1.25e-03 1.14e-02 | ord  2.17  2.02 -3.19
rw_t10       sqv  be       X f1     lev 32-64-128-256-512  L1 2.01e-01 1.48e-01 1.37e-01 7.04e-02 | ord  0.44  0.11  0.96
rw_t10       sqv  be       X dens   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 4.32e-04 | ord  2.18  2.14  1.62
rw_t10       sqv  be       X velx   lev 32-64-128-256-512  L1 2.72e-02 6.07e-03 1.37e-03 4.65e-04 | ord  2.17  2.15  1.56
rw_t10       sqv  be       X eint   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 4.52e-04 | ord  2.18  2.14  1.56
rw_t10       sqv  be       T e      lev 1-2-4-8-16  L1 7.83e-02 4.28e-02 2.24e-02 1.15e-02 | ord  0.87  0.93  0.96
rw_t10       sqv  be       T f1     lev 1-2-4-8-16  L1 6.68e-02 3.57e-02 1.85e-02 9.37e-03 | ord  0.90  0.95  0.98
rw_t10       sqv  be       T dens   lev 1-2-4-8-16  L1 6.76e-03 3.37e-03 1.68e-03 8.41e-04 | ord  1.00  1.00  1.00
rw_t10       sqv  be       T velx   lev 1-2-4-8-16  L1 6.41e-03 3.20e-03 1.60e-03 8.01e-04 | ord  1.00  1.00  1.00
rw_t10       sqv  be       T eint   lev 1-2-4-8-16  L1 6.79e-03 3.39e-03 1.69e-03 8.45e-04 | ord  1.00  1.00  1.00
rw_t10       sqv  be       C e      lev 32-64-128-256-512  L1 2.03e-01 1.38e-01 7.98e-02 5.34e-02 | ord  0.56  0.79  0.58
rw_t10       sqv  be       C f1     lev 32-64-128-256-512  L1 2.09e-01 1.24e-01 6.77e-02 6.49e-02 | ord  0.76  0.87  0.06
rw_t10       sqv  be       C dens   lev 32-64-128-256-512  L1 4.39e-02 1.52e-02 6.87e-03 3.53e-03 | ord  1.53  1.15  0.96
rw_t10       sqv  be       C velx   lev 32-64-128-256-512  L1 4.01e-02 1.37e-02 6.48e-03 3.33e-03 | ord  1.55  1.08  0.96
rw_t10       sqv  be       C eint   lev 32-64-128-256-512  L1 4.44e-02 1.55e-02 6.93e-03 3.56e-03 | ord  1.52  1.16  0.96
rw_t10       sqv  hesdirk2 X e      lev 32-64-128-256-512  L1 1.54e-02 4.10e-03 1.05e-03 2.85e-04 | ord  1.91  1.97  1.88
rw_t10       sqv  hesdirk2 X f1     lev 32-64-128-256-512  L1 1.68e-01 4.76e-02 4.51e-02 2.29e-03 | ord  1.82  0.08  4.30
rw_t10       sqv  hesdirk2 X dens   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.93e-04 | ord  2.18  2.13  2.21
rw_t10       sqv  hesdirk2 X velx   lev 32-64-128-256-512  L1 2.74e-02 6.13e-03 1.39e-03 3.12e-04 | ord  2.16  2.14  2.16
rw_t10       sqv  hesdirk2 X eint   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.93e-04 | ord  2.18  2.13  2.21
rw_t10       sqv  hesdirk2 T e      lev 1-2-4-8-16  L1 8.55e-04 4.27e-04 4.75e-04 5.20e-04 | ord  1.00 -0.15 -0.13
rw_t10       sqv  hesdirk2 T f1     lev 1-2-4-8-16  L1 4.40e-04 4.94e-04 2.33e-04 8.51e-05 | ord -0.17  1.08  1.45
rw_t10       sqv  hesdirk2 T dens   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       sqv  hesdirk2 T velx   lev 1-2-4-8-16  L1 1.29e-03 6.37e-04 3.17e-04 1.58e-04 | ord  1.02  1.01  1.01
rw_t10       sqv  hesdirk2 T eint   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       sqv  hesdirk2 C e      lev 32-64-128-256-512  L1 1.43e-02 2.81e-03 5.28e-04 3.16e-04 | ord  2.35  2.41  0.74
rw_t10       sqv  hesdirk2 C f1     lev 32-64-128-256-512  L1 1.57e-01 2.82e-02 4.02e-03 2.24e-03 | ord  2.48  2.81  0.84
rw_t10       sqv  hesdirk2 C dens   lev 32-64-128-256-512  L1 2.86e-02 6.61e-03 1.70e-03 6.13e-04 | ord  2.12  1.96  1.47
rw_t10       sqv  hesdirk2 C velx   lev 32-64-128-256-512  L1 2.87e-02 6.43e-03 1.53e-03 6.16e-04 | ord  2.16  2.07  1.31
rw_t10       sqv  hesdirk2 C eint   lev 32-64-128-256-512  L1 2.87e-02 6.61e-03 1.70e-03 6.13e-04 | ord  2.12  1.96  1.47
rw_t10       sqa  be       X e      lev 32-64-128-256-512  L1 2.28e-02 5.07e-03 1.25e-03 1.14e-02 | ord  2.17  2.02 -3.19
rw_t10       sqa  be       X f1     lev 32-64-128-256-512  L1 2.01e-01 1.48e-01 1.37e-01 7.04e-02 | ord  0.44  0.11  0.96
rw_t10       sqa  be       X dens   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 4.32e-04 | ord  2.18  2.14  1.62
rw_t10       sqa  be       X velx   lev 32-64-128-256-512  L1 2.72e-02 6.07e-03 1.37e-03 4.65e-04 | ord  2.17  2.15  1.56
rw_t10       sqa  be       X eint   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 4.52e-04 | ord  2.18  2.14  1.56
rw_t10       sqa  be       T e      lev 1-2-4-8-16  L1 7.83e-02 4.28e-02 2.24e-02 1.15e-02 | ord  0.87  0.93  0.96
rw_t10       sqa  be       T f1     lev 1-2-4-8-16  L1 6.68e-02 3.57e-02 1.85e-02 9.37e-03 | ord  0.90  0.95  0.98
rw_t10       sqa  be       T dens   lev 1-2-4-8-16  L1 6.76e-03 3.37e-03 1.68e-03 8.41e-04 | ord  1.00  1.00  1.00
rw_t10       sqa  be       T velx   lev 1-2-4-8-16  L1 6.41e-03 3.20e-03 1.60e-03 8.01e-04 | ord  1.00  1.00  1.00
rw_t10       sqa  be       T eint   lev 1-2-4-8-16  L1 6.79e-03 3.39e-03 1.69e-03 8.45e-04 | ord  1.00  1.00  1.00
rw_t10       sqa  be       C e      lev 32-64-128-256-512  L1 2.03e-01 1.38e-01 7.98e-02 5.34e-02 | ord  0.56  0.79  0.58
rw_t10       sqa  be       C f1     lev 32-64-128-256-512  L1 2.09e-01 1.24e-01 6.77e-02 6.49e-02 | ord  0.76  0.87  0.06
rw_t10       sqa  be       C dens   lev 32-64-128-256-512  L1 4.39e-02 1.52e-02 6.87e-03 3.53e-03 | ord  1.53  1.15  0.96
rw_t10       sqa  be       C velx   lev 32-64-128-256-512  L1 4.01e-02 1.37e-02 6.48e-03 3.33e-03 | ord  1.55  1.08  0.96
rw_t10       sqa  be       C eint   lev 32-64-128-256-512  L1 4.44e-02 1.55e-02 6.93e-03 3.56e-03 | ord  1.52  1.16  0.96
rw_t10       sqa  hesdirk2 X e      lev 32-64-128-256-512  L1 1.54e-02 4.10e-03 1.05e-03 2.85e-04 | ord  1.91  1.97  1.88
rw_t10       sqa  hesdirk2 X f1     lev 32-64-128-256-512  L1 1.68e-01 4.76e-02 4.51e-02 2.29e-03 | ord  1.82  0.08  4.30
rw_t10       sqa  hesdirk2 X dens   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.93e-04 | ord  2.18  2.13  2.21
rw_t10       sqa  hesdirk2 X velx   lev 32-64-128-256-512  L1 2.74e-02 6.13e-03 1.39e-03 3.12e-04 | ord  2.16  2.14  2.16
rw_t10       sqa  hesdirk2 X eint   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.93e-04 | ord  2.18  2.13  2.21
rw_t10       sqa  hesdirk2 T e      lev 1-2-4-8-16  L1 8.55e-04 4.27e-04 4.75e-04 5.20e-04 | ord  1.00 -0.15 -0.13
rw_t10       sqa  hesdirk2 T f1     lev 1-2-4-8-16  L1 4.40e-04 4.94e-04 2.33e-04 8.51e-05 | ord -0.17  1.08  1.45
rw_t10       sqa  hesdirk2 T dens   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       sqa  hesdirk2 T velx   lev 1-2-4-8-16  L1 1.29e-03 6.37e-04 3.17e-04 1.58e-04 | ord  1.02  1.01  1.01
rw_t10       sqa  hesdirk2 T eint   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       sqa  hesdirk2 C e      lev 32-64-128-256-512  L1 1.43e-02 2.81e-03 5.28e-04 3.16e-04 | ord  2.35  2.41  0.74
rw_t10       sqa  hesdirk2 C f1     lev 32-64-128-256-512  L1 1.57e-01 2.82e-02 4.02e-03 2.24e-03 | ord  2.48  2.81  0.84
rw_t10       sqa  hesdirk2 C dens   lev 32-64-128-256-512  L1 2.86e-02 6.61e-03 1.70e-03 6.13e-04 | ord  2.12  1.96  1.47
rw_t10       sqa  hesdirk2 C velx   lev 32-64-128-256-512  L1 2.87e-02 6.43e-03 1.53e-03 6.16e-04 | ord  2.16  2.07  1.31
rw_t10       sqa  hesdirk2 C eint   lev 32-64-128-256-512  L1 2.87e-02 6.61e-03 1.70e-03 6.13e-04 | ord  2.12  1.96  1.47
rw_t10       aq   be       X e      lev 32-64-128-256-512  L1 2.26e-02 3.61e-02 1.10e-01 1.96e-01 | ord -0.68 -1.61 -0.83
rw_t10       aq   be       X f1     lev 32-64-128-256-512  L1 2.05e-01 2.38e-01 5.72e-01 1.01e+00 | ord -0.21 -1.26 -0.82
rw_t10       aq   be       X dens   lev 32-64-128-256-512  L1 2.66e-02 6.21e-03 2.64e-03 2.76e-03 | ord  2.10  1.23 -0.06
rw_t10       aq   be       X velx   lev 32-64-128-256-512  L1 2.73e-02 6.41e-03 3.36e-03 4.84e-03 | ord  2.09  0.93 -0.53
rw_t10       aq   be       X eint   lev 32-64-128-256-512  L1 2.66e-02 6.26e-03 2.82e-03 3.02e-03 | ord  2.09  1.15 -0.10
rw_t10       aq   be       T e      lev 1-2-4-8-16  L1 8.35e-02 3.47e-02 1.19e-02 1.23e-03 | ord  1.27  1.54  3.28
rw_t10       aq   be       T f1     lev 1-2-4-8-16  L1 6.71e-02 4.99e-02 6.60e-02 7.30e-02 | ord  0.43 -0.40 -0.15
rw_t10       aq   be       T dens   lev 1-2-4-8-16  L1 5.03e-03 2.55e-03 1.31e-03 7.01e-04 | ord  0.98  0.95  0.91
rw_t10       aq   be       T velx   lev 1-2-4-8-16  L1 5.22e-03 2.64e-03 1.36e-03 7.46e-04 | ord  0.98  0.95  0.87
rw_t10       aq   be       T eint   lev 1-2-4-8-16  L1 5.05e-03 2.55e-03 1.32e-03 7.01e-04 | ord  0.98  0.96  0.91
rw_t10       aq   be       C e      lev 32-64-128-256-512  L1 2.61e-01 2.60e-01 2.09e-01 3.24e-02 | ord  0.01  0.31  2.69
rw_t10       aq   be       C f1     lev 32-64-128-256-512  L1 2.56e-01 2.83e-01 7.11e-01 1.05e+00 | ord -0.14 -1.33 -0.57
rw_t10       aq   be       C dens   lev 32-64-128-256-512  L1 4.52e-02 1.71e-02 8.07e-03 2.60e-03 | ord  1.40  1.08  1.64
rw_t10       aq   be       C velx   lev 32-64-128-256-512  L1 4.12e-02 1.51e-02 7.71e-03 4.32e-03 | ord  1.44  0.97  0.83
rw_t10       aq   be       C eint   lev 32-64-128-256-512  L1 4.58e-02 1.75e-02 8.22e-03 2.60e-03 | ord  1.39  1.09  1.66
rw_t10       aq   hesdirk2 X e      lev 32-64-128-256-512  L1 1.56e-02 7.20e-03 2.00e-02 4.00e-02 | ord  1.11 -1.47 -1.00
rw_t10       aq   hesdirk2 X f1     lev 32-64-128-256-512  L1 1.68e-01 1.52e-01 4.40e-01 8.61e-01 | ord  0.14 -1.53 -0.97
rw_t10       aq   hesdirk2 X dens   lev 32-64-128-256-512  L1 2.67e-02 5.85e-03 1.20e-03 5.23e-04 | ord  2.19  2.29  1.20
rw_t10       aq   hesdirk2 X velx   lev 32-64-128-256-512  L1 2.74e-02 6.10e-03 2.38e-03 4.59e-03 | ord  2.17  1.36 -0.95
rw_t10       aq   hesdirk2 X eint   lev 32-64-128-256-512  L1 2.67e-02 5.84e-03 1.17e-03 5.51e-04 | ord  2.19  2.32  1.09
rw_t10       aq   hesdirk2 T e      lev 1-2-4-8-16  L1 2.06e-03 2.90e-03 2.87e-03 2.22e-03 | ord -0.49  0.02  0.37
rw_t10       aq   hesdirk2 T f1     lev 1-2-4-8-16  L1 6.11e-02 7.23e-02 6.80e-02 5.17e-02 | ord -0.24  0.09  0.40
rw_t10       aq   hesdirk2 T dens   lev 1-2-4-8-16  L1 1.30e-03 6.66e-04 3.48e-04 1.82e-04 | ord  0.96  0.94  0.93
rw_t10       aq   hesdirk2 T velx   lev 1-2-4-8-16  L1 1.32e-03 7.05e-04 4.19e-04 2.71e-04 | ord  0.91  0.75  0.63
rw_t10       aq   hesdirk2 T eint   lev 1-2-4-8-16  L1 1.30e-03 6.67e-04 3.48e-04 1.83e-04 | ord  0.96  0.94  0.93
rw_t10       aq   hesdirk2 C e      lev 32-64-128-256-512  L1 1.57e-02 1.17e-02 2.28e-02 3.17e-02 | ord  0.42 -0.96 -0.47
rw_t10       aq   hesdirk2 C f1     lev 32-64-128-256-512  L1 1.78e-01 2.27e-01 4.98e-01 6.71e-01 | ord -0.35 -1.13 -0.43
rw_t10       aq   hesdirk2 C dens   lev 32-64-128-256-512  L1 2.86e-02 6.49e-03 1.39e-03 2.78e-04 | ord  2.14  2.22  2.33
rw_t10       aq   hesdirk2 C velx   lev 32-64-128-256-512  L1 2.86e-02 6.38e-03 2.77e-03 3.61e-03 | ord  2.17  1.20 -0.38
rw_t10       aq   hesdirk2 C eint   lev 32-64-128-256-512  L1 2.86e-02 6.48e-03 1.36e-03 3.08e-04 | ord  2.14  2.25  2.15
rw_t10       kn01 be       X e      lev 32-64-128-256-512  L1 2.28e-02 5.07e-03 1.25e-03 3.83e-04 | ord  2.17  2.02  1.70
rw_t10       kn01 be       X f1     lev 32-64-128-256-512  L1 2.01e-01 1.48e-01 1.37e-01 3.33e-03 | ord  0.44  0.12  5.36
rw_t10       kn01 be       X dens   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 2.98e-04 | ord  2.18  2.14  2.16
rw_t10       kn01 be       X velx   lev 32-64-128-256-512  L1 2.72e-02 6.07e-03 1.37e-03 3.13e-04 | ord  2.17  2.15  2.13
rw_t10       kn01 be       X eint   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 2.98e-04 | ord  2.18  2.14  2.16
rw_t10       kn01 be       T e      lev 1-2-4-8-16  L1 7.84e-02 4.28e-02 2.24e-02 1.15e-02 | ord  0.87  0.93  0.96
rw_t10       kn01 be       T f1     lev 1-2-4-8-16  L1 6.67e-02 3.57e-02 1.86e-02 9.29e-03 | ord  0.90  0.94  1.00
rw_t10       kn01 be       T dens   lev 1-2-4-8-16  L1 6.76e-03 3.37e-03 1.68e-03 8.41e-04 | ord  1.00  1.00  1.00
rw_t10       kn01 be       T velx   lev 1-2-4-8-16  L1 6.41e-03 3.20e-03 1.60e-03 8.01e-04 | ord  1.00  1.00  1.00
rw_t10       kn01 be       T eint   lev 1-2-4-8-16  L1 6.79e-03 3.39e-03 1.69e-03 8.45e-04 | ord  1.00  1.00  1.00
rw_t10       kn01 be       C e      lev 32-64-128-256-512  L1 2.04e-01 1.38e-01 7.95e-02 4.36e-02 | ord  0.56  0.79  0.87
rw_t10       kn01 be       C f1     lev 32-64-128-256-512  L1 2.09e-01 1.24e-01 6.74e-02 3.64e-02 | ord  0.76  0.88  0.89
rw_t10       kn01 be       C dens   lev 32-64-128-256-512  L1 4.39e-02 1.52e-02 6.86e-03 3.40e-03 | ord  1.53  1.15  1.01
rw_t10       kn01 be       C velx   lev 32-64-128-256-512  L1 4.01e-02 1.37e-02 6.47e-03 3.22e-03 | ord  1.55  1.08  1.01
rw_t10       kn01 be       C eint   lev 32-64-128-256-512  L1 4.44e-02 1.55e-02 6.92e-03 3.42e-03 | ord  1.52  1.16  1.02
rw_t10       kn01 hesdirk2 X e      lev 32-64-128-256-512  L1 1.54e-02 4.10e-03 1.05e-03 3.24e-04 | ord  1.91  1.97  1.69
rw_t10       kn01 hesdirk2 X f1     lev 32-64-128-256-512  L1 1.68e-01 4.76e-02 4.51e-02 1.75e-03 | ord  1.82  0.08  4.69
rw_t10       kn01 hesdirk2 X dens   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.93e-04 | ord  2.18  2.13  2.20
rw_t10       kn01 hesdirk2 X velx   lev 32-64-128-256-512  L1 2.74e-02 6.13e-03 1.39e-03 3.13e-04 | ord  2.16  2.14  2.15
rw_t10       kn01 hesdirk2 X eint   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.94e-04 | ord  2.18  2.13  2.20
rw_t10       kn01 hesdirk2 T e      lev 1-2-4-8-16  L1 8.55e-04 4.29e-04 2.16e-04 1.06e-04 | ord  1.00  0.99  1.03
rw_t10       kn01 hesdirk2 T f1     lev 1-2-4-8-16  L1 4.62e-04 5.83e-04 2.01e-04 1.23e-04 | ord -0.33  1.54  0.70
rw_t10       kn01 hesdirk2 T dens   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       kn01 hesdirk2 T velx   lev 1-2-4-8-16  L1 1.29e-03 6.37e-04 3.17e-04 1.58e-04 | ord  1.02  1.01  1.01
rw_t10       kn01 hesdirk2 T eint   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       kn01 hesdirk2 C e      lev 32-64-128-256-512  L1 1.43e-02 2.81e-03 5.27e-04 2.49e-04 | ord  2.35  2.41  1.08
rw_t10       kn01 hesdirk2 C f1     lev 32-64-128-256-512  L1 1.57e-01 2.81e-02 3.77e-03 1.90e-03 | ord  2.48  2.90  0.99
rw_t10       kn01 hesdirk2 C dens   lev 32-64-128-256-512  L1 2.86e-02 6.61e-03 1.70e-03 6.14e-04 | ord  2.12  1.96  1.47
rw_t10       kn01 hesdirk2 C velx   lev 32-64-128-256-512  L1 2.87e-02 6.43e-03 1.53e-03 6.17e-04 | ord  2.16  2.07  1.31
rw_t10       kn01 hesdirk2 C eint   lev 32-64-128-256-512  L1 2.87e-02 6.61e-03 1.70e-03 6.14e-04 | ord  2.12  1.96  1.47
rw_t10       kn03 be       X e      lev 32-64-128-256-512  L1 2.28e-02 5.07e-03 1.25e-03 3.83e-04 | ord  2.17  2.02  1.70
rw_t10       kn03 be       X f1     lev 32-64-128-256-512  L1 2.01e-01 1.48e-01 1.37e-01 3.33e-03 | ord  0.44  0.12  5.36
rw_t10       kn03 be       X dens   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 2.98e-04 | ord  2.18  2.14  2.16
rw_t10       kn03 be       X velx   lev 32-64-128-256-512  L1 2.72e-02 6.07e-03 1.37e-03 3.13e-04 | ord  2.17  2.15  2.13
rw_t10       kn03 be       X eint   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 2.98e-04 | ord  2.18  2.14  2.16
rw_t10       kn03 be       T e      lev 1-2-4-8-16  L1 7.84e-02 4.28e-02 2.24e-02 1.15e-02 | ord  0.87  0.93  0.96
rw_t10       kn03 be       T f1     lev 1-2-4-8-16  L1 6.67e-02 3.57e-02 1.86e-02 9.29e-03 | ord  0.90  0.94  1.00
rw_t10       kn03 be       T dens   lev 1-2-4-8-16  L1 6.76e-03 3.37e-03 1.68e-03 8.41e-04 | ord  1.00  1.00  1.00
rw_t10       kn03 be       T velx   lev 1-2-4-8-16  L1 6.41e-03 3.20e-03 1.60e-03 8.01e-04 | ord  1.00  1.00  1.00
rw_t10       kn03 be       T eint   lev 1-2-4-8-16  L1 6.79e-03 3.39e-03 1.69e-03 8.45e-04 | ord  1.00  1.00  1.00
rw_t10       kn03 be       C e      lev 32-64-128-256-512  L1 2.04e-01 1.38e-01 7.95e-02 4.36e-02 | ord  0.56  0.79  0.87
rw_t10       kn03 be       C f1     lev 32-64-128-256-512  L1 2.09e-01 1.24e-01 6.74e-02 3.64e-02 | ord  0.76  0.88  0.89
rw_t10       kn03 be       C dens   lev 32-64-128-256-512  L1 4.39e-02 1.52e-02 6.86e-03 3.40e-03 | ord  1.53  1.15  1.01
rw_t10       kn03 be       C velx   lev 32-64-128-256-512  L1 4.01e-02 1.37e-02 6.47e-03 3.22e-03 | ord  1.55  1.08  1.01
rw_t10       kn03 be       C eint   lev 32-64-128-256-512  L1 4.44e-02 1.55e-02 6.92e-03 3.42e-03 | ord  1.52  1.16  1.02
rw_t10       kn03 hesdirk2 X e      lev 32-64-128-256-512  L1 1.54e-02 4.10e-03 1.05e-03 3.24e-04 | ord  1.91  1.97  1.69
rw_t10       kn03 hesdirk2 X f1     lev 32-64-128-256-512  L1 1.68e-01 4.76e-02 4.51e-02 1.75e-03 | ord  1.82  0.08  4.69
rw_t10       kn03 hesdirk2 X dens   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.93e-04 | ord  2.18  2.13  2.20
rw_t10       kn03 hesdirk2 X velx   lev 32-64-128-256-512  L1 2.74e-02 6.13e-03 1.39e-03 3.13e-04 | ord  2.16  2.14  2.15
rw_t10       kn03 hesdirk2 X eint   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.94e-04 | ord  2.18  2.13  2.20
rw_t10       kn03 hesdirk2 T e      lev 1-2-4-8-16  L1 8.55e-04 4.29e-04 2.16e-04 1.06e-04 | ord  1.00  0.99  1.03
rw_t10       kn03 hesdirk2 T f1     lev 1-2-4-8-16  L1 4.62e-04 5.83e-04 2.01e-04 1.23e-04 | ord -0.33  1.54  0.70
rw_t10       kn03 hesdirk2 T dens   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       kn03 hesdirk2 T velx   lev 1-2-4-8-16  L1 1.29e-03 6.37e-04 3.17e-04 1.58e-04 | ord  1.02  1.01  1.01
rw_t10       kn03 hesdirk2 T eint   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       kn03 hesdirk2 C e      lev 32-64-128-256-512  L1 1.43e-02 2.81e-03 5.27e-04 2.49e-04 | ord  2.35  2.41  1.08
rw_t10       kn03 hesdirk2 C f1     lev 32-64-128-256-512  L1 1.57e-01 2.81e-02 3.77e-03 1.90e-03 | ord  2.48  2.90  0.99
rw_t10       kn03 hesdirk2 C dens   lev 32-64-128-256-512  L1 2.86e-02 6.61e-03 1.70e-03 6.14e-04 | ord  2.12  1.96  1.47
rw_t10       kn03 hesdirk2 C velx   lev 32-64-128-256-512  L1 2.87e-02 6.43e-03 1.53e-03 6.17e-04 | ord  2.16  2.07  1.31
rw_t10       kn03 hesdirk2 C eint   lev 32-64-128-256-512  L1 2.87e-02 6.61e-03 1.70e-03 6.14e-04 | ord  2.12  1.96  1.47
rw_t10       kn1  be       X e      lev 32-64-128-256-512  L1 2.28e-02 5.07e-03 1.25e-03 3.83e-04 | ord  2.17  2.02  1.70
rw_t10       kn1  be       X f1     lev 32-64-128-256-512  L1 2.01e-01 1.48e-01 1.37e-01 3.33e-03 | ord  0.44  0.12  5.36
rw_t10       kn1  be       X dens   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 2.98e-04 | ord  2.18  2.14  2.16
rw_t10       kn1  be       X velx   lev 32-64-128-256-512  L1 2.72e-02 6.07e-03 1.37e-03 3.13e-04 | ord  2.17  2.15  2.13
rw_t10       kn1  be       X eint   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 2.98e-04 | ord  2.18  2.14  2.16
rw_t10       kn1  be       T e      lev 1-2-4-8-16  L1 7.84e-02 4.28e-02 2.24e-02 1.15e-02 | ord  0.87  0.93  0.96
rw_t10       kn1  be       T f1     lev 1-2-4-8-16  L1 6.67e-02 3.57e-02 1.86e-02 9.29e-03 | ord  0.90  0.94  1.00
rw_t10       kn1  be       T dens   lev 1-2-4-8-16  L1 6.76e-03 3.37e-03 1.68e-03 8.41e-04 | ord  1.00  1.00  1.00
rw_t10       kn1  be       T velx   lev 1-2-4-8-16  L1 6.41e-03 3.20e-03 1.60e-03 8.01e-04 | ord  1.00  1.00  1.00
rw_t10       kn1  be       T eint   lev 1-2-4-8-16  L1 6.79e-03 3.39e-03 1.69e-03 8.45e-04 | ord  1.00  1.00  1.00
rw_t10       kn1  be       C e      lev 32-64-128-256-512  L1 2.04e-01 1.38e-01 7.95e-02 4.36e-02 | ord  0.56  0.79  0.87
rw_t10       kn1  be       C f1     lev 32-64-128-256-512  L1 2.09e-01 1.24e-01 6.74e-02 3.64e-02 | ord  0.76  0.88  0.89
rw_t10       kn1  be       C dens   lev 32-64-128-256-512  L1 4.39e-02 1.52e-02 6.86e-03 3.40e-03 | ord  1.53  1.15  1.01
rw_t10       kn1  be       C velx   lev 32-64-128-256-512  L1 4.01e-02 1.37e-02 6.47e-03 3.22e-03 | ord  1.55  1.08  1.01
rw_t10       kn1  be       C eint   lev 32-64-128-256-512  L1 4.44e-02 1.55e-02 6.92e-03 3.42e-03 | ord  1.52  1.16  1.02
rw_t10       kn1  hesdirk2 X e      lev 32-64-128-256-512  L1 1.54e-02 4.10e-03 1.05e-03 3.24e-04 | ord  1.91  1.97  1.69
rw_t10       kn1  hesdirk2 X f1     lev 32-64-128-256-512  L1 1.68e-01 4.76e-02 4.51e-02 1.75e-03 | ord  1.82  0.08  4.69
rw_t10       kn1  hesdirk2 X dens   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.93e-04 | ord  2.18  2.13  2.20
rw_t10       kn1  hesdirk2 X velx   lev 32-64-128-256-512  L1 2.74e-02 6.13e-03 1.39e-03 3.13e-04 | ord  2.16  2.14  2.15
rw_t10       kn1  hesdirk2 X eint   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.94e-04 | ord  2.18  2.13  2.20
rw_t10       kn1  hesdirk2 T e      lev 1-2-4-8-16  L1 8.55e-04 4.29e-04 2.16e-04 1.06e-04 | ord  1.00  0.99  1.03
rw_t10       kn1  hesdirk2 T f1     lev 1-2-4-8-16  L1 4.62e-04 5.83e-04 2.01e-04 1.23e-04 | ord -0.33  1.54  0.70
rw_t10       kn1  hesdirk2 T dens   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       kn1  hesdirk2 T velx   lev 1-2-4-8-16  L1 1.29e-03 6.37e-04 3.17e-04 1.58e-04 | ord  1.02  1.01  1.01
rw_t10       kn1  hesdirk2 T eint   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       kn1  hesdirk2 C e      lev 32-64-128-256-512  L1 1.43e-02 2.81e-03 5.27e-04 2.49e-04 | ord  2.35  2.41  1.08
rw_t10       kn1  hesdirk2 C f1     lev 32-64-128-256-512  L1 1.57e-01 2.81e-02 3.77e-03 1.90e-03 | ord  2.48  2.90  0.99
rw_t10       kn1  hesdirk2 C dens   lev 32-64-128-256-512  L1 2.86e-02 6.61e-03 1.70e-03 6.14e-04 | ord  2.12  1.96  1.47
rw_t10       kn1  hesdirk2 C velx   lev 32-64-128-256-512  L1 2.87e-02 6.43e-03 1.53e-03 6.17e-04 | ord  2.16  2.07  1.31
rw_t10       kn1  hesdirk2 C eint   lev 32-64-128-256-512  L1 2.87e-02 6.61e-03 1.70e-03 6.14e-04 | ord  2.12  1.96  1.47
rw_t10       knp  be       X e      lev 32-64-128-256-512  L1 2.28e-02 5.07e-03 1.25e-03 3.83e-04 | ord  2.17  2.02  1.70
rw_t10       knp  be       X f1     lev 32-64-128-256-512  L1 2.01e-01 1.48e-01 1.37e-01 3.33e-03 | ord  0.44  0.12  5.36
rw_t10       knp  be       X dens   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 2.98e-04 | ord  2.18  2.14  2.16
rw_t10       knp  be       X velx   lev 32-64-128-256-512  L1 2.72e-02 6.07e-03 1.37e-03 3.13e-04 | ord  2.17  2.15  2.13
rw_t10       knp  be       X eint   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 2.98e-04 | ord  2.18  2.14  2.16
rw_t10       knp  be       T e      lev 1-2-4-8-16  L1 7.84e-02 4.28e-02 2.24e-02 1.15e-02 | ord  0.87  0.93  0.96
rw_t10       knp  be       T f1     lev 1-2-4-8-16  L1 6.67e-02 3.57e-02 1.86e-02 9.29e-03 | ord  0.90  0.94  1.00
rw_t10       knp  be       T dens   lev 1-2-4-8-16  L1 6.76e-03 3.37e-03 1.68e-03 8.41e-04 | ord  1.00  1.00  1.00
rw_t10       knp  be       T velx   lev 1-2-4-8-16  L1 6.41e-03 3.20e-03 1.60e-03 8.01e-04 | ord  1.00  1.00  1.00
rw_t10       knp  be       T eint   lev 1-2-4-8-16  L1 6.79e-03 3.39e-03 1.69e-03 8.45e-04 | ord  1.00  1.00  1.00
rw_t10       knp  be       C e      lev 32-64-128-256-512  L1 2.04e-01 1.38e-01 7.95e-02 4.36e-02 | ord  0.56  0.79  0.87
rw_t10       knp  be       C f1     lev 32-64-128-256-512  L1 2.09e-01 1.24e-01 6.74e-02 3.64e-02 | ord  0.76  0.88  0.89
rw_t10       knp  be       C dens   lev 32-64-128-256-512  L1 4.39e-02 1.52e-02 6.86e-03 3.40e-03 | ord  1.53  1.15  1.01
rw_t10       knp  be       C velx   lev 32-64-128-256-512  L1 4.01e-02 1.37e-02 6.47e-03 3.22e-03 | ord  1.55  1.08  1.01
rw_t10       knp  be       C eint   lev 32-64-128-256-512  L1 4.44e-02 1.55e-02 6.92e-03 3.42e-03 | ord  1.52  1.16  1.02
rw_t10       knp  hesdirk2 X e      lev 32-64-128-256-512  L1 1.54e-02 4.10e-03 1.05e-03 3.24e-04 | ord  1.91  1.97  1.69
rw_t10       knp  hesdirk2 X f1     lev 32-64-128-256-512  L1 1.68e-01 4.76e-02 4.51e-02 1.75e-03 | ord  1.82  0.08  4.69
rw_t10       knp  hesdirk2 X dens   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.93e-04 | ord  2.18  2.13  2.20
rw_t10       knp  hesdirk2 X velx   lev 32-64-128-256-512  L1 2.74e-02 6.13e-03 1.39e-03 3.13e-04 | ord  2.16  2.14  2.15
rw_t10       knp  hesdirk2 X eint   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.94e-04 | ord  2.18  2.13  2.20
rw_t10       knp  hesdirk2 T e      lev 1-2-4-8-16  L1 8.55e-04 4.29e-04 2.16e-04 1.06e-04 | ord  1.00  0.99  1.03
rw_t10       knp  hesdirk2 T f1     lev 1-2-4-8-16  L1 4.62e-04 5.83e-04 2.01e-04 1.23e-04 | ord -0.33  1.54  0.70
rw_t10       knp  hesdirk2 T dens   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       knp  hesdirk2 T velx   lev 1-2-4-8-16  L1 1.29e-03 6.37e-04 3.17e-04 1.58e-04 | ord  1.02  1.01  1.01
rw_t10       knp  hesdirk2 T eint   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       knp  hesdirk2 C e      lev 32-64-128-256-512  L1 1.43e-02 2.81e-03 5.27e-04 2.49e-04 | ord  2.35  2.41  1.08
rw_t10       knp  hesdirk2 C f1     lev 32-64-128-256-512  L1 1.57e-01 2.81e-02 3.77e-03 1.90e-03 | ord  2.48  2.90  0.99
rw_t10       knp  hesdirk2 C dens   lev 32-64-128-256-512  L1 2.86e-02 6.61e-03 1.70e-03 6.14e-04 | ord  2.12  1.96  1.47
rw_t10       knp  hesdirk2 C velx   lev 32-64-128-256-512  L1 2.87e-02 6.43e-03 1.53e-03 6.17e-04 | ord  2.16  2.07  1.31
rw_t10       knp  hesdirk2 C eint   lev 32-64-128-256-512  L1 2.87e-02 6.61e-03 1.70e-03 6.14e-04 | ord  2.12  1.96  1.47
rw_t1000     cen  be       X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     cen  be       X f1     lev 32-64-128-256  L1 4.93e-03 1.47e-02 5.91e-04 | ord -1.58  4.64
rw_t1000     cen  be       X dens   lev 32-64-128-256  L1 8.73e-03 2.04e-03 5.31e-04 | ord  2.10  1.94
rw_t1000     cen  be       X velx   lev 32-64-128-256  L1 8.40e-03 2.03e-03 5.36e-04 | ord  2.05  1.92
rw_t1000     cen  be       X eint   lev 32-64-128-256  L1 8.69e-03 2.02e-03 5.30e-04 | ord  2.10  1.93
rw_t1000     cen  be       T e      lev 2-4-8-16  L1 3.30e-03 1.68e-03 8.50e-04 | ord  0.97  0.99
rw_t1000     cen  be       T f1     lev 2-4-8-16  L1 2.77e-02 1.69e-03 8.52e-04 | ord  4.03  0.99
rw_t1000     cen  be       T dens   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.41e-04 | ord  0.98  0.99
rw_t1000     cen  be       T velx   lev 2-4-8-16  L1 3.31e-03 1.68e-03 8.50e-04 | ord  0.97  0.99
rw_t1000     cen  be       T eint   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.43e-04 | ord  0.97  0.99
rw_t1000     cen  be       C e      lev 32-64-128-256  L1 2.49e-02 1.15e-02 5.58e-03 | ord  1.12  1.04
rw_t1000     cen  be       C f1     lev 32-64-128-256  L1 2.62e-02 1.31e-02 2.99e-02 | ord  1.00 -1.19
rw_t1000     cen  be       C dens   lev 32-64-128-256  L1 2.62e-02 1.17e-02 5.62e-03 | ord  1.16  1.06
rw_t1000     cen  be       C velx   lev 32-64-128-256  L1 2.51e-02 1.17e-02 5.67e-03 | ord  1.10  1.05
rw_t1000     cen  be       C eint   lev 32-64-128-256  L1 2.56e-02 1.16e-02 5.60e-03 | ord  1.14  1.05
rw_t1000     cen  hesdirk2 X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     cen  hesdirk2 X f1     lev 32-64-128-256  L1 5.75e-03 5.87e-03 5.91e-04 | ord -0.03  3.31
rw_t1000     cen  hesdirk2 X dens   lev 32-64-128-256  L1 8.76e-03 2.04e-03 5.27e-04 | ord  2.11  1.95
rw_t1000     cen  hesdirk2 X velx   lev 32-64-128-256  L1 8.44e-03 2.04e-03 5.37e-04 | ord  2.05  1.92
rw_t1000     cen  hesdirk2 X eint   lev 32-64-128-256  L1 8.72e-03 2.01e-03 5.27e-04 | ord  2.12  1.93
rw_t1000     cen  hesdirk2 T e      lev 2-4-8-16  L1 1.14e-02 5.66e-03 2.81e-03 | ord  1.01  1.01
rw_t1000     cen  hesdirk2 T f1     lev 2-4-8-16  L1 1.20e-02 5.54e-03 2.79e-03 | ord  1.12  0.99
rw_t1000     cen  hesdirk2 T dens   lev 2-4-8-16  L1 1.16e-02 5.78e-03 2.88e-03 | ord  1.01  1.00
rw_t1000     cen  hesdirk2 T velx   lev 2-4-8-16  L1 1.09e-02 5.54e-03 2.78e-03 | ord  0.98  0.99
rw_t1000     cen  hesdirk2 T eint   lev 2-4-8-16  L1 1.14e-02 5.68e-03 2.83e-03 | ord  1.01  1.00
rw_t1000     cen  hesdirk2 C e      lev 32-64-128-256  L1 7.99e-02 3.81e-02 1.84e-02 | ord  1.07  1.05
rw_t1000     cen  hesdirk2 C f1     lev 32-64-128-256  L1 6.32e-02 3.35e-02 1.78e-02 | ord  0.92  0.91
rw_t1000     cen  hesdirk2 C dens   lev 32-64-128-256  L1 7.95e-02 3.85e-02 1.88e-02 | ord  1.05  1.04
rw_t1000     cen  hesdirk2 C velx   lev 32-64-128-256  L1 6.43e-02 3.36e-02 1.72e-02 | ord  0.94  0.96
rw_t1000     cen  hesdirk2 C eint   lev 32-64-128-256  L1 7.92e-02 3.82e-02 1.85e-02 | ord  1.05  1.04
rw_t1000     hr   be       X e      lev 32-64-128-256  L1 8.56e-03 2.20e-03 8.67e-04 | ord  1.96  1.34
rw_t1000     hr   be       X f1     lev 32-64-128-256  L1 4.89e-03 1.49e-02 9.42e-04 | ord -1.61  3.98
rw_t1000     hr   be       X dens   lev 32-64-128-256  L1 8.81e-03 2.20e-03 8.60e-04 | ord  2.00  1.36
rw_t1000     hr   be       X velx   lev 32-64-128-256  L1 8.48e-03 2.20e-03 8.68e-04 | ord  1.95  1.34
rw_t1000     hr   be       X eint   lev 32-64-128-256  L1 8.78e-03 2.19e-03 8.59e-04 | ord  2.00  1.35
rw_t1000     hr   be       T e      lev 2-4-8-16  L1 3.30e-03 1.68e-03 8.46e-04 | ord  0.97  0.99
rw_t1000     hr   be       T f1     lev 2-4-8-16  L1 2.77e-02 1.69e-03 8.48e-04 | ord  4.03  0.99
rw_t1000     hr   be       T dens   lev 2-4-8-16  L1 3.27e-03 1.66e-03 8.37e-04 | ord  0.98  0.99
rw_t1000     hr   be       T velx   lev 2-4-8-16  L1 3.31e-03 1.68e-03 8.47e-04 | ord  0.97  0.99
rw_t1000     hr   be       T eint   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.40e-04 | ord  0.97  0.99
rw_t1000     hr   be       C e      lev 32-64-128-256  L1 2.49e-02 1.16e-02 5.91e-03 | ord  1.10  0.98
rw_t1000     hr   be       C f1     lev 32-64-128-256  L1 2.63e-02 1.33e-02 3.03e-02 | ord  0.98 -1.19
rw_t1000     hr   be       C dens   lev 32-64-128-256  L1 2.63e-02 1.19e-02 5.95e-03 | ord  1.14  1.00
rw_t1000     hr   be       C velx   lev 32-64-128-256  L1 2.52e-02 1.19e-02 6.00e-03 | ord  1.08  0.99
rw_t1000     hr   be       C eint   lev 32-64-128-256  L1 2.57e-02 1.18e-02 5.94e-03 | ord  1.12  0.99
rw_t1000     hr   hesdirk2 X e      lev 32-64-128-256  L1 8.56e-03 2.20e-03 8.62e-04 | ord  1.96  1.35
rw_t1000     hr   hesdirk2 X f1     lev 32-64-128-256  L1 5.84e-03 6.04e-03 9.37e-04 | ord -0.05  2.69
rw_t1000     hr   hesdirk2 X dens   lev 32-64-128-256  L1 8.84e-03 2.20e-03 8.51e-04 | ord  2.01  1.37
rw_t1000     hr   hesdirk2 X velx   lev 32-64-128-256  L1 8.52e-03 2.20e-03 8.64e-04 | ord  1.95  1.35
rw_t1000     hr   hesdirk2 X eint   lev 32-64-128-256  L1 8.80e-03 2.17e-03 8.52e-04 | ord  2.02  1.35
rw_t1000     hr   hesdirk2 T e      lev 2-4-8-16  L1 1.14e-02 5.66e-03 2.81e-03 | ord  1.01  1.01
rw_t1000     hr   hesdirk2 T f1     lev 2-4-8-16  L1 1.20e-02 5.54e-03 2.78e-03 | ord  1.12  0.99
rw_t1000     hr   hesdirk2 T dens   lev 2-4-8-16  L1 1.16e-02 5.78e-03 2.88e-03 | ord  1.01  1.00
rw_t1000     hr   hesdirk2 T velx   lev 2-4-8-16  L1 1.09e-02 5.54e-03 2.78e-03 | ord  0.98  0.99
rw_t1000     hr   hesdirk2 T eint   lev 2-4-8-16  L1 1.14e-02 5.68e-03 2.83e-03 | ord  1.01  1.01
rw_t1000     hr   hesdirk2 C e      lev 32-64-128-256  L1 8.00e-02 3.81e-02 1.84e-02 | ord  1.07  1.05
rw_t1000     hr   hesdirk2 C f1     lev 32-64-128-256  L1 6.33e-02 3.35e-02 1.79e-02 | ord  0.92  0.90
rw_t1000     hr   hesdirk2 C dens   lev 32-64-128-256  L1 7.95e-02 3.86e-02 1.89e-02 | ord  1.04  1.03
rw_t1000     hr   hesdirk2 C velx   lev 32-64-128-256  L1 6.43e-02 3.36e-02 1.73e-02 | ord  0.94  0.96
rw_t1000     hr   hesdirk2 C eint   lev 32-64-128-256  L1 7.93e-02 3.82e-02 1.86e-02 | ord  1.05  1.04
rw_t1000     sq   be       X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     sq   be       X f1     lev 32-64-128-256  L1 4.93e-03 1.47e-02 5.91e-04 | ord -1.58  4.64
rw_t1000     sq   be       X dens   lev 32-64-128-256  L1 8.73e-03 2.04e-03 5.31e-04 | ord  2.10  1.94
rw_t1000     sq   be       X velx   lev 32-64-128-256  L1 8.40e-03 2.03e-03 5.36e-04 | ord  2.05  1.92
rw_t1000     sq   be       X eint   lev 32-64-128-256  L1 8.69e-03 2.02e-03 5.30e-04 | ord  2.10  1.93
rw_t1000     sq   be       T e      lev 2-4-8-16  L1 3.30e-03 1.68e-03 8.50e-04 | ord  0.97  0.98
rw_t1000     sq   be       T f1     lev 2-4-8-16  L1 2.77e-02 1.69e-03 8.52e-04 | ord  4.03  0.99
rw_t1000     sq   be       T dens   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.41e-04 | ord  0.98  0.99
rw_t1000     sq   be       T velx   lev 2-4-8-16  L1 3.31e-03 1.68e-03 8.50e-04 | ord  0.97  0.99
rw_t1000     sq   be       T eint   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.43e-04 | ord  0.97  0.99
rw_t1000     sq   be       C e      lev 32-64-128-256  L1 2.49e-02 1.15e-02 5.58e-03 | ord  1.12  1.04
rw_t1000     sq   be       C f1     lev 32-64-128-256  L1 2.62e-02 1.31e-02 2.99e-02 | ord  1.00 -1.19
rw_t1000     sq   be       C dens   lev 32-64-128-256  L1 2.62e-02 1.17e-02 5.62e-03 | ord  1.16  1.06
rw_t1000     sq   be       C velx   lev 32-64-128-256  L1 2.51e-02 1.17e-02 5.67e-03 | ord  1.10  1.05
rw_t1000     sq   be       C eint   lev 32-64-128-256  L1 2.56e-02 1.16e-02 5.60e-03 | ord  1.14  1.05
rw_t1000     sq   hesdirk2 X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     sq   hesdirk2 X f1     lev 32-64-128-256  L1 5.75e-03 5.87e-03 5.91e-04 | ord -0.03  3.31
rw_t1000     sq   hesdirk2 X dens   lev 32-64-128-256  L1 8.76e-03 2.04e-03 5.27e-04 | ord  2.11  1.95
rw_t1000     sq   hesdirk2 X velx   lev 32-64-128-256  L1 8.44e-03 2.04e-03 5.37e-04 | ord  2.05  1.92
rw_t1000     sq   hesdirk2 X eint   lev 32-64-128-256  L1 8.72e-03 2.01e-03 5.27e-04 | ord  2.12  1.93
rw_t1000     sq   hesdirk2 T e      lev 2-4-8-16  L1 1.14e-02 5.66e-03 2.81e-03 | ord  1.01  1.01
rw_t1000     sq   hesdirk2 T f1     lev 2-4-8-16  L1 1.20e-02 5.54e-03 2.79e-03 | ord  1.12  0.99
rw_t1000     sq   hesdirk2 T dens   lev 2-4-8-16  L1 1.16e-02 5.78e-03 2.88e-03 | ord  1.01  1.00
rw_t1000     sq   hesdirk2 T velx   lev 2-4-8-16  L1 1.09e-02 5.54e-03 2.78e-03 | ord  0.98  0.99
rw_t1000     sq   hesdirk2 T eint   lev 2-4-8-16  L1 1.14e-02 5.68e-03 2.83e-03 | ord  1.01  1.00
rw_t1000     sq   hesdirk2 C e      lev 32-64-128-256  L1 7.99e-02 3.81e-02 1.84e-02 | ord  1.07  1.05
rw_t1000     sq   hesdirk2 C f1     lev 32-64-128-256  L1 6.32e-02 3.35e-02 1.78e-02 | ord  0.92  0.91
rw_t1000     sq   hesdirk2 C dens   lev 32-64-128-256  L1 7.95e-02 3.85e-02 1.88e-02 | ord  1.05  1.04
rw_t1000     sq   hesdirk2 C velx   lev 32-64-128-256  L1 6.43e-02 3.36e-02 1.72e-02 | ord  0.94  0.96
rw_t1000     sq   hesdirk2 C eint   lev 32-64-128-256  L1 7.92e-02 3.82e-02 1.85e-02 | ord  1.05  1.04
rw_t1000     sqf  be       X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     sqf  be       X f1     lev 32-64-128-256  L1 4.93e-03 1.47e-02 5.91e-04 | ord -1.58  4.64
rw_t1000     sqf  be       X dens   lev 32-64-128-256  L1 8.73e-03 2.04e-03 5.31e-04 | ord  2.10  1.94
rw_t1000     sqf  be       X velx   lev 32-64-128-256  L1 8.40e-03 2.03e-03 5.36e-04 | ord  2.05  1.92
rw_t1000     sqf  be       X eint   lev 32-64-128-256  L1 8.69e-03 2.02e-03 5.30e-04 | ord  2.10  1.93
rw_t1000     sqf  be       T e      lev 2-4-8-16  L1 3.30e-03 1.68e-03 8.50e-04 | ord  0.97  0.98
rw_t1000     sqf  be       T f1     lev 2-4-8-16  L1 2.77e-02 1.69e-03 8.52e-04 | ord  4.03  0.99
rw_t1000     sqf  be       T dens   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.41e-04 | ord  0.98  0.99
rw_t1000     sqf  be       T velx   lev 2-4-8-16  L1 3.31e-03 1.68e-03 8.50e-04 | ord  0.97  0.99
rw_t1000     sqf  be       T eint   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.43e-04 | ord  0.97  0.99
rw_t1000     sqf  be       C e      lev 32-64-128-256  L1 2.49e-02 1.15e-02 5.58e-03 | ord  1.12  1.04
rw_t1000     sqf  be       C f1     lev 32-64-128-256  L1 2.62e-02 1.31e-02 2.99e-02 | ord  1.00 -1.19
rw_t1000     sqf  be       C dens   lev 32-64-128-256  L1 2.62e-02 1.17e-02 5.62e-03 | ord  1.16  1.06
rw_t1000     sqf  be       C velx   lev 32-64-128-256  L1 2.51e-02 1.17e-02 5.67e-03 | ord  1.10  1.05
rw_t1000     sqf  be       C eint   lev 32-64-128-256  L1 2.56e-02 1.16e-02 5.60e-03 | ord  1.14  1.05
rw_t1000     sqf  hesdirk2 X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     sqf  hesdirk2 X f1     lev 32-64-128-256  L1 5.75e-03 5.87e-03 5.91e-04 | ord -0.03  3.31
rw_t1000     sqf  hesdirk2 X dens   lev 32-64-128-256  L1 8.76e-03 2.04e-03 5.27e-04 | ord  2.11  1.95
rw_t1000     sqf  hesdirk2 X velx   lev 32-64-128-256  L1 8.44e-03 2.04e-03 5.37e-04 | ord  2.05  1.92
rw_t1000     sqf  hesdirk2 X eint   lev 32-64-128-256  L1 8.72e-03 2.01e-03 5.27e-04 | ord  2.12  1.93
rw_t1000     sqf  hesdirk2 T e      lev 2-4-8-16  L1 1.14e-02 5.66e-03 2.81e-03 | ord  1.01  1.01
rw_t1000     sqf  hesdirk2 T f1     lev 2-4-8-16  L1 1.20e-02 5.54e-03 2.79e-03 | ord  1.12  0.99
rw_t1000     sqf  hesdirk2 T dens   lev 2-4-8-16  L1 1.16e-02 5.78e-03 2.88e-03 | ord  1.01  1.00
rw_t1000     sqf  hesdirk2 T velx   lev 2-4-8-16  L1 1.09e-02 5.54e-03 2.78e-03 | ord  0.98  0.99
rw_t1000     sqf  hesdirk2 T eint   lev 2-4-8-16  L1 1.14e-02 5.68e-03 2.83e-03 | ord  1.01  1.00
rw_t1000     sqf  hesdirk2 C e      lev 32-64-128-256  L1 7.99e-02 3.81e-02 1.84e-02 | ord  1.07  1.05
rw_t1000     sqf  hesdirk2 C f1     lev 32-64-128-256  L1 6.32e-02 3.35e-02 1.78e-02 | ord  0.92  0.91
rw_t1000     sqf  hesdirk2 C dens   lev 32-64-128-256  L1 7.95e-02 3.85e-02 1.88e-02 | ord  1.05  1.04
rw_t1000     sqf  hesdirk2 C velx   lev 32-64-128-256  L1 6.43e-02 3.36e-02 1.72e-02 | ord  0.94  0.96
rw_t1000     sqf  hesdirk2 C eint   lev 32-64-128-256  L1 7.92e-02 3.82e-02 1.85e-02 | ord  1.05  1.04
rw_t1000     sqv  be       X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     sqv  be       X f1     lev 32-64-128-256  L1 4.93e-03 1.47e-02 5.91e-04 | ord -1.58  4.64
rw_t1000     sqv  be       X dens   lev 32-64-128-256  L1 8.73e-03 2.04e-03 5.31e-04 | ord  2.10  1.94
rw_t1000     sqv  be       X velx   lev 32-64-128-256  L1 8.40e-03 2.03e-03 5.36e-04 | ord  2.05  1.92
rw_t1000     sqv  be       X eint   lev 32-64-128-256  L1 8.69e-03 2.02e-03 5.30e-04 | ord  2.10  1.93
rw_t1000     sqv  be       T e      lev 2-4-8-16  L1 3.30e-03 1.68e-03 8.50e-04 | ord  0.97  0.98
rw_t1000     sqv  be       T f1     lev 2-4-8-16  L1 2.77e-02 1.69e-03 8.52e-04 | ord  4.03  0.99
rw_t1000     sqv  be       T dens   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.41e-04 | ord  0.98  0.99
rw_t1000     sqv  be       T velx   lev 2-4-8-16  L1 3.31e-03 1.68e-03 8.50e-04 | ord  0.97  0.99
rw_t1000     sqv  be       T eint   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.43e-04 | ord  0.97  0.99
rw_t1000     sqv  be       C e      lev 32-64-128-256  L1 2.49e-02 1.15e-02 5.58e-03 | ord  1.12  1.04
rw_t1000     sqv  be       C f1     lev 32-64-128-256  L1 2.62e-02 1.31e-02 2.99e-02 | ord  1.00 -1.19
rw_t1000     sqv  be       C dens   lev 32-64-128-256  L1 2.62e-02 1.17e-02 5.62e-03 | ord  1.16  1.06
rw_t1000     sqv  be       C velx   lev 32-64-128-256  L1 2.51e-02 1.17e-02 5.67e-03 | ord  1.10  1.05
rw_t1000     sqv  be       C eint   lev 32-64-128-256  L1 2.56e-02 1.16e-02 5.60e-03 | ord  1.14  1.05
rw_t1000     sqv  hesdirk2 X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     sqv  hesdirk2 X f1     lev 32-64-128-256  L1 5.75e-03 5.87e-03 5.91e-04 | ord -0.03  3.31
rw_t1000     sqv  hesdirk2 X dens   lev 32-64-128-256  L1 8.76e-03 2.04e-03 5.27e-04 | ord  2.11  1.95
rw_t1000     sqv  hesdirk2 X velx   lev 32-64-128-256  L1 8.44e-03 2.04e-03 5.37e-04 | ord  2.05  1.92
rw_t1000     sqv  hesdirk2 X eint   lev 32-64-128-256  L1 8.72e-03 2.01e-03 5.27e-04 | ord  2.12  1.93
rw_t1000     sqv  hesdirk2 T e      lev 2-4-8-16  L1 1.14e-02 5.66e-03 2.81e-03 | ord  1.01  1.01
rw_t1000     sqv  hesdirk2 T f1     lev 2-4-8-16  L1 1.20e-02 5.54e-03 2.79e-03 | ord  1.12  0.99
rw_t1000     sqv  hesdirk2 T dens   lev 2-4-8-16  L1 1.16e-02 5.78e-03 2.88e-03 | ord  1.01  1.00
rw_t1000     sqv  hesdirk2 T velx   lev 2-4-8-16  L1 1.09e-02 5.54e-03 2.78e-03 | ord  0.98  0.99
rw_t1000     sqv  hesdirk2 T eint   lev 2-4-8-16  L1 1.14e-02 5.68e-03 2.83e-03 | ord  1.01  1.00
rw_t1000     sqv  hesdirk2 C e      lev 32-64-128-256  L1 7.99e-02 3.81e-02 1.84e-02 | ord  1.07  1.05
rw_t1000     sqv  hesdirk2 C f1     lev 32-64-128-256  L1 6.32e-02 3.35e-02 1.78e-02 | ord  0.92  0.91
rw_t1000     sqv  hesdirk2 C dens   lev 32-64-128-256  L1 7.95e-02 3.85e-02 1.88e-02 | ord  1.05  1.04
rw_t1000     sqv  hesdirk2 C velx   lev 32-64-128-256  L1 6.43e-02 3.36e-02 1.72e-02 | ord  0.94  0.96
rw_t1000     sqv  hesdirk2 C eint   lev 32-64-128-256  L1 7.92e-02 3.82e-02 1.85e-02 | ord  1.05  1.04
rw_t1000     sqa  be       X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     sqa  be       X f1     lev 32-64-128-256  L1 4.93e-03 1.47e-02 5.91e-04 | ord -1.58  4.64
rw_t1000     sqa  be       X dens   lev 32-64-128-256  L1 8.73e-03 2.04e-03 5.31e-04 | ord  2.10  1.94
rw_t1000     sqa  be       X velx   lev 32-64-128-256  L1 8.40e-03 2.03e-03 5.36e-04 | ord  2.05  1.92
rw_t1000     sqa  be       X eint   lev 32-64-128-256  L1 8.69e-03 2.02e-03 5.30e-04 | ord  2.10  1.93
rw_t1000     sqa  be       T e      lev 2-4-8-16  L1 3.30e-03 1.68e-03 8.50e-04 | ord  0.97  0.98
rw_t1000     sqa  be       T f1     lev 2-4-8-16  L1 2.77e-02 1.69e-03 8.52e-04 | ord  4.03  0.99
rw_t1000     sqa  be       T dens   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.41e-04 | ord  0.98  0.99
rw_t1000     sqa  be       T velx   lev 2-4-8-16  L1 3.31e-03 1.68e-03 8.50e-04 | ord  0.97  0.99
rw_t1000     sqa  be       T eint   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.43e-04 | ord  0.97  0.99
rw_t1000     sqa  be       C e      lev 32-64-128-256  L1 2.49e-02 1.15e-02 5.58e-03 | ord  1.12  1.04
rw_t1000     sqa  be       C f1     lev 32-64-128-256  L1 2.62e-02 1.31e-02 2.99e-02 | ord  1.00 -1.19
rw_t1000     sqa  be       C dens   lev 32-64-128-256  L1 2.62e-02 1.17e-02 5.62e-03 | ord  1.16  1.06
rw_t1000     sqa  be       C velx   lev 32-64-128-256  L1 2.51e-02 1.17e-02 5.67e-03 | ord  1.10  1.05
rw_t1000     sqa  be       C eint   lev 32-64-128-256  L1 2.56e-02 1.16e-02 5.60e-03 | ord  1.14  1.05
rw_t1000     sqa  hesdirk2 X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     sqa  hesdirk2 X f1     lev 32-64-128-256  L1 5.75e-03 5.87e-03 5.91e-04 | ord -0.03  3.31
rw_t1000     sqa  hesdirk2 X dens   lev 32-64-128-256  L1 8.76e-03 2.04e-03 5.27e-04 | ord  2.11  1.95
rw_t1000     sqa  hesdirk2 X velx   lev 32-64-128-256  L1 8.44e-03 2.04e-03 5.37e-04 | ord  2.05  1.92
rw_t1000     sqa  hesdirk2 X eint   lev 32-64-128-256  L1 8.72e-03 2.01e-03 5.27e-04 | ord  2.12  1.93
rw_t1000     sqa  hesdirk2 T e      lev 2-4-8-16  L1 1.14e-02 5.66e-03 2.81e-03 | ord  1.01  1.01
rw_t1000     sqa  hesdirk2 T f1     lev 2-4-8-16  L1 1.20e-02 5.54e-03 2.79e-03 | ord  1.12  0.99
rw_t1000     sqa  hesdirk2 T dens   lev 2-4-8-16  L1 1.16e-02 5.78e-03 2.88e-03 | ord  1.01  1.00
rw_t1000     sqa  hesdirk2 T velx   lev 2-4-8-16  L1 1.09e-02 5.54e-03 2.78e-03 | ord  0.98  0.99
rw_t1000     sqa  hesdirk2 T eint   lev 2-4-8-16  L1 1.14e-02 5.68e-03 2.83e-03 | ord  1.01  1.00
rw_t1000     sqa  hesdirk2 C e      lev 32-64-128-256  L1 7.99e-02 3.81e-02 1.84e-02 | ord  1.07  1.05
rw_t1000     sqa  hesdirk2 C f1     lev 32-64-128-256  L1 6.32e-02 3.35e-02 1.78e-02 | ord  0.92  0.91
rw_t1000     sqa  hesdirk2 C dens   lev 32-64-128-256  L1 7.95e-02 3.85e-02 1.88e-02 | ord  1.05  1.04
rw_t1000     sqa  hesdirk2 C velx   lev 32-64-128-256  L1 6.43e-02 3.36e-02 1.72e-02 | ord  0.94  0.96
rw_t1000     sqa  hesdirk2 C eint   lev 32-64-128-256  L1 7.92e-02 3.82e-02 1.85e-02 | ord  1.05  1.04
rw_t1000     aq   be       X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.37e-04 | ord  2.06  1.92
rw_t1000     aq   be       X f1     lev 32-64-128-256  L1 4.93e-03 1.47e-02 5.92e-04 | ord -1.58  4.63
rw_t1000     aq   be       X dens   lev 32-64-128-256  L1 8.72e-03 2.04e-03 5.32e-04 | ord  2.10  1.94
rw_t1000     aq   be       X velx   lev 32-64-128-256  L1 8.40e-03 2.03e-03 5.37e-04 | ord  2.05  1.92
rw_t1000     aq   be       X eint   lev 32-64-128-256  L1 8.69e-03 2.02e-03 5.31e-04 | ord  2.10  1.93
rw_t1000     aq   be       T e      lev 2-4-8-16  L1 3.30e-03 1.68e-03 8.50e-04 | ord  0.97  0.98
rw_t1000     aq   be       T f1     lev 2-4-8-16  L1 2.77e-02 1.69e-03 8.52e-04 | ord  4.03  0.99
rw_t1000     aq   be       T dens   lev 2-4-8-16  L1 3.27e-03 1.67e-03 8.41e-04 | ord  0.98  0.99
rw_t1000     aq   be       T velx   lev 2-4-8-16  L1 3.31e-03 1.68e-03 8.50e-04 | ord  0.97  0.99
rw_t1000     aq   be       T eint   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.43e-04 | ord  0.97  0.99
rw_t1000     aq   be       C e      lev 32-64-128-256  L1 2.49e-02 1.15e-02 5.58e-03 | ord  1.12  1.04
rw_t1000     aq   be       C f1     lev 32-64-128-256  L1 2.62e-02 1.31e-02 2.99e-02 | ord  1.00 -1.19
rw_t1000     aq   be       C dens   lev 32-64-128-256  L1 2.62e-02 1.17e-02 5.62e-03 | ord  1.16  1.06
rw_t1000     aq   be       C velx   lev 32-64-128-256  L1 2.51e-02 1.17e-02 5.67e-03 | ord  1.10  1.05
rw_t1000     aq   be       C eint   lev 32-64-128-256  L1 2.56e-02 1.16e-02 5.60e-03 | ord  1.14  1.05
rw_t1000     aq   hesdirk2 X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     aq   hesdirk2 X f1     lev 32-64-128-256  L1 5.75e-03 5.87e-03 5.91e-04 | ord -0.03  3.31
rw_t1000     aq   hesdirk2 X dens   lev 32-64-128-256  L1 8.76e-03 2.04e-03 5.27e-04 | ord  2.11  1.95
rw_t1000     aq   hesdirk2 X velx   lev 32-64-128-256  L1 8.44e-03 2.04e-03 5.37e-04 | ord  2.05  1.92
rw_t1000     aq   hesdirk2 X eint   lev 32-64-128-256  L1 8.72e-03 2.01e-03 5.26e-04 | ord  2.12  1.93
rw_t1000     aq   hesdirk2 T e      lev 2-4-8-16  L1 1.14e-02 5.66e-03 2.81e-03 | ord  1.01  1.01
rw_t1000     aq   hesdirk2 T f1     lev 2-4-8-16  L1 1.20e-02 5.54e-03 2.79e-03 | ord  1.12  0.99
rw_t1000     aq   hesdirk2 T dens   lev 2-4-8-16  L1 1.16e-02 5.78e-03 2.88e-03 | ord  1.01  1.00
rw_t1000     aq   hesdirk2 T velx   lev 2-4-8-16  L1 1.09e-02 5.54e-03 2.79e-03 | ord  0.98  0.99
rw_t1000     aq   hesdirk2 T eint   lev 2-4-8-16  L1 1.14e-02 5.68e-03 2.83e-03 | ord  1.01  1.00
rw_t1000     aq   hesdirk2 C e      lev 32-64-128-256  L1 7.99e-02 3.81e-02 1.84e-02 | ord  1.07  1.05
rw_t1000     aq   hesdirk2 C f1     lev 32-64-128-256  L1 6.32e-02 3.35e-02 1.78e-02 | ord  0.92  0.92
rw_t1000     aq   hesdirk2 C dens   lev 32-64-128-256  L1 7.95e-02 3.85e-02 1.88e-02 | ord  1.05  1.04
rw_t1000     aq   hesdirk2 C velx   lev 32-64-128-256  L1 6.43e-02 3.36e-02 1.72e-02 | ord  0.94  0.96
rw_t1000     aq   hesdirk2 C eint   lev 32-64-128-256  L1 7.92e-02 3.82e-02 1.85e-02 | ord  1.05  1.04
rw_t1000     kn01 be       X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     kn01 be       X f1     lev 32-64-128-256  L1 4.93e-03 1.47e-02 5.91e-04 | ord -1.58  4.64
rw_t1000     kn01 be       X dens   lev 32-64-128-256  L1 8.73e-03 2.04e-03 5.31e-04 | ord  2.10  1.94
rw_t1000     kn01 be       X velx   lev 32-64-128-256  L1 8.40e-03 2.03e-03 5.36e-04 | ord  2.05  1.92
rw_t1000     kn01 be       X eint   lev 32-64-128-256  L1 8.69e-03 2.02e-03 5.30e-04 | ord  2.10  1.93
rw_t1000     kn01 be       T e      lev 2-4-8-16  L1 3.30e-03 1.68e-03 8.50e-04 | ord  0.97  0.98
rw_t1000     kn01 be       T f1     lev 2-4-8-16  L1 2.77e-02 1.69e-03 8.52e-04 | ord  4.03  0.99
rw_t1000     kn01 be       T dens   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.41e-04 | ord  0.98  0.99
rw_t1000     kn01 be       T velx   lev 2-4-8-16  L1 3.31e-03 1.68e-03 8.50e-04 | ord  0.97  0.99
rw_t1000     kn01 be       T eint   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.43e-04 | ord  0.97  0.99
rw_t1000     kn01 be       C e      lev 32-64-128-256  L1 2.49e-02 1.15e-02 5.58e-03 | ord  1.12  1.04
rw_t1000     kn01 be       C f1     lev 32-64-128-256  L1 2.62e-02 1.31e-02 2.99e-02 | ord  1.00 -1.19
rw_t1000     kn01 be       C dens   lev 32-64-128-256  L1 2.62e-02 1.17e-02 5.62e-03 | ord  1.16  1.06
rw_t1000     kn01 be       C velx   lev 32-64-128-256  L1 2.51e-02 1.17e-02 5.67e-03 | ord  1.10  1.05
rw_t1000     kn01 be       C eint   lev 32-64-128-256  L1 2.56e-02 1.16e-02 5.60e-03 | ord  1.14  1.05
rw_t1000     kn01 hesdirk2 X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     kn01 hesdirk2 X f1     lev 32-64-128-256  L1 5.75e-03 5.87e-03 5.91e-04 | ord -0.03  3.31
rw_t1000     kn01 hesdirk2 X dens   lev 32-64-128-256  L1 8.76e-03 2.04e-03 5.27e-04 | ord  2.11  1.95
rw_t1000     kn01 hesdirk2 X velx   lev 32-64-128-256  L1 8.44e-03 2.04e-03 5.37e-04 | ord  2.05  1.92
rw_t1000     kn01 hesdirk2 X eint   lev 32-64-128-256  L1 8.72e-03 2.01e-03 5.27e-04 | ord  2.12  1.93
rw_t1000     kn01 hesdirk2 T e      lev 2-4-8-16  L1 1.14e-02 5.66e-03 2.81e-03 | ord  1.01  1.01
rw_t1000     kn01 hesdirk2 T f1     lev 2-4-8-16  L1 1.20e-02 5.54e-03 2.79e-03 | ord  1.12  0.99
rw_t1000     kn01 hesdirk2 T dens   lev 2-4-8-16  L1 1.16e-02 5.78e-03 2.88e-03 | ord  1.01  1.00
rw_t1000     kn01 hesdirk2 T velx   lev 2-4-8-16  L1 1.09e-02 5.54e-03 2.78e-03 | ord  0.98  0.99
rw_t1000     kn01 hesdirk2 T eint   lev 2-4-8-16  L1 1.14e-02 5.68e-03 2.83e-03 | ord  1.01  1.00
rw_t1000     kn01 hesdirk2 C e      lev 32-64-128-256  L1 7.99e-02 3.81e-02 1.84e-02 | ord  1.07  1.05
rw_t1000     kn01 hesdirk2 C f1     lev 32-64-128-256  L1 6.32e-02 3.35e-02 1.78e-02 | ord  0.92  0.91
rw_t1000     kn01 hesdirk2 C dens   lev 32-64-128-256  L1 7.95e-02 3.85e-02 1.88e-02 | ord  1.05  1.04
rw_t1000     kn01 hesdirk2 C velx   lev 32-64-128-256  L1 6.43e-02 3.36e-02 1.72e-02 | ord  0.94  0.96
rw_t1000     kn01 hesdirk2 C eint   lev 32-64-128-256  L1 7.92e-02 3.82e-02 1.85e-02 | ord  1.05  1.04
rw_t1000     kn03 be       X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     kn03 be       X f1     lev 32-64-128-256  L1 4.93e-03 1.47e-02 5.91e-04 | ord -1.58  4.64
rw_t1000     kn03 be       X dens   lev 32-64-128-256  L1 8.73e-03 2.04e-03 5.31e-04 | ord  2.10  1.94
rw_t1000     kn03 be       X velx   lev 32-64-128-256  L1 8.40e-03 2.03e-03 5.36e-04 | ord  2.05  1.92
rw_t1000     kn03 be       X eint   lev 32-64-128-256  L1 8.69e-03 2.02e-03 5.30e-04 | ord  2.10  1.93
rw_t1000     kn03 be       T e      lev 2-4-8-16  L1 3.30e-03 1.68e-03 8.50e-04 | ord  0.97  0.98
rw_t1000     kn03 be       T f1     lev 2-4-8-16  L1 2.77e-02 1.69e-03 8.52e-04 | ord  4.03  0.99
rw_t1000     kn03 be       T dens   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.41e-04 | ord  0.98  0.99
rw_t1000     kn03 be       T velx   lev 2-4-8-16  L1 3.31e-03 1.68e-03 8.50e-04 | ord  0.97  0.99
rw_t1000     kn03 be       T eint   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.43e-04 | ord  0.97  0.99
rw_t1000     kn03 be       C e      lev 32-64-128-256  L1 2.49e-02 1.15e-02 5.58e-03 | ord  1.12  1.04
rw_t1000     kn03 be       C f1     lev 32-64-128-256  L1 2.62e-02 1.31e-02 2.99e-02 | ord  1.00 -1.19
rw_t1000     kn03 be       C dens   lev 32-64-128-256  L1 2.62e-02 1.17e-02 5.62e-03 | ord  1.16  1.06
rw_t1000     kn03 be       C velx   lev 32-64-128-256  L1 2.51e-02 1.17e-02 5.67e-03 | ord  1.10  1.05
rw_t1000     kn03 be       C eint   lev 32-64-128-256  L1 2.56e-02 1.16e-02 5.60e-03 | ord  1.14  1.05
rw_t1000     kn03 hesdirk2 X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     kn03 hesdirk2 X f1     lev 32-64-128-256  L1 5.75e-03 5.87e-03 5.91e-04 | ord -0.03  3.31
rw_t1000     kn03 hesdirk2 X dens   lev 32-64-128-256  L1 8.76e-03 2.04e-03 5.27e-04 | ord  2.11  1.95
rw_t1000     kn03 hesdirk2 X velx   lev 32-64-128-256  L1 8.44e-03 2.04e-03 5.37e-04 | ord  2.05  1.92
rw_t1000     kn03 hesdirk2 X eint   lev 32-64-128-256  L1 8.72e-03 2.01e-03 5.27e-04 | ord  2.12  1.93
rw_t1000     kn03 hesdirk2 T e      lev 2-4-8-16  L1 1.14e-02 5.66e-03 2.81e-03 | ord  1.01  1.01
rw_t1000     kn03 hesdirk2 T f1     lev 2-4-8-16  L1 1.20e-02 5.54e-03 2.79e-03 | ord  1.12  0.99
rw_t1000     kn03 hesdirk2 T dens   lev 2-4-8-16  L1 1.16e-02 5.78e-03 2.88e-03 | ord  1.01  1.00
rw_t1000     kn03 hesdirk2 T velx   lev 2-4-8-16  L1 1.09e-02 5.54e-03 2.78e-03 | ord  0.98  0.99
rw_t1000     kn03 hesdirk2 T eint   lev 2-4-8-16  L1 1.14e-02 5.68e-03 2.83e-03 | ord  1.01  1.00
rw_t1000     kn03 hesdirk2 C e      lev 32-64-128-256  L1 7.99e-02 3.81e-02 1.84e-02 | ord  1.07  1.05
rw_t1000     kn03 hesdirk2 C f1     lev 32-64-128-256  L1 6.32e-02 3.35e-02 1.78e-02 | ord  0.92  0.91
rw_t1000     kn03 hesdirk2 C dens   lev 32-64-128-256  L1 7.95e-02 3.85e-02 1.88e-02 | ord  1.05  1.04
rw_t1000     kn03 hesdirk2 C velx   lev 32-64-128-256  L1 6.43e-02 3.36e-02 1.72e-02 | ord  0.94  0.96
rw_t1000     kn03 hesdirk2 C eint   lev 32-64-128-256  L1 7.92e-02 3.82e-02 1.85e-02 | ord  1.05  1.04
rw_t1000     kn1  be       X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     kn1  be       X f1     lev 32-64-128-256  L1 4.93e-03 1.47e-02 5.91e-04 | ord -1.58  4.64
rw_t1000     kn1  be       X dens   lev 32-64-128-256  L1 8.73e-03 2.04e-03 5.31e-04 | ord  2.10  1.94
rw_t1000     kn1  be       X velx   lev 32-64-128-256  L1 8.40e-03 2.03e-03 5.36e-04 | ord  2.05  1.92
rw_t1000     kn1  be       X eint   lev 32-64-128-256  L1 8.69e-03 2.02e-03 5.30e-04 | ord  2.10  1.93
rw_t1000     kn1  be       T e      lev 2-4-8-16  L1 3.30e-03 1.68e-03 8.50e-04 | ord  0.97  0.98
rw_t1000     kn1  be       T f1     lev 2-4-8-16  L1 2.77e-02 1.69e-03 8.52e-04 | ord  4.03  0.99
rw_t1000     kn1  be       T dens   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.41e-04 | ord  0.98  0.99
rw_t1000     kn1  be       T velx   lev 2-4-8-16  L1 3.31e-03 1.68e-03 8.50e-04 | ord  0.97  0.99
rw_t1000     kn1  be       T eint   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.43e-04 | ord  0.97  0.99
rw_t1000     kn1  be       C e      lev 32-64-128-256  L1 2.49e-02 1.15e-02 5.58e-03 | ord  1.12  1.04
rw_t1000     kn1  be       C f1     lev 32-64-128-256  L1 2.62e-02 1.31e-02 2.99e-02 | ord  1.00 -1.19
rw_t1000     kn1  be       C dens   lev 32-64-128-256  L1 2.62e-02 1.17e-02 5.62e-03 | ord  1.16  1.06
rw_t1000     kn1  be       C velx   lev 32-64-128-256  L1 2.51e-02 1.17e-02 5.67e-03 | ord  1.10  1.05
rw_t1000     kn1  be       C eint   lev 32-64-128-256  L1 2.56e-02 1.16e-02 5.60e-03 | ord  1.14  1.05
rw_t1000     kn1  hesdirk2 X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     kn1  hesdirk2 X f1     lev 32-64-128-256  L1 5.75e-03 5.87e-03 5.91e-04 | ord -0.03  3.31
rw_t1000     kn1  hesdirk2 X dens   lev 32-64-128-256  L1 8.76e-03 2.04e-03 5.27e-04 | ord  2.11  1.95
rw_t1000     kn1  hesdirk2 X velx   lev 32-64-128-256  L1 8.44e-03 2.04e-03 5.37e-04 | ord  2.05  1.92
rw_t1000     kn1  hesdirk2 X eint   lev 32-64-128-256  L1 8.72e-03 2.01e-03 5.27e-04 | ord  2.12  1.93
rw_t1000     kn1  hesdirk2 T e      lev 2-4-8-16  L1 1.14e-02 5.66e-03 2.81e-03 | ord  1.01  1.01
rw_t1000     kn1  hesdirk2 T f1     lev 2-4-8-16  L1 1.20e-02 5.54e-03 2.79e-03 | ord  1.12  0.99
rw_t1000     kn1  hesdirk2 T dens   lev 2-4-8-16  L1 1.16e-02 5.78e-03 2.88e-03 | ord  1.01  1.00
rw_t1000     kn1  hesdirk2 T velx   lev 2-4-8-16  L1 1.09e-02 5.54e-03 2.78e-03 | ord  0.98  0.99
rw_t1000     kn1  hesdirk2 T eint   lev 2-4-8-16  L1 1.14e-02 5.68e-03 2.83e-03 | ord  1.01  1.00
rw_t1000     kn1  hesdirk2 C e      lev 32-64-128-256  L1 7.99e-02 3.81e-02 1.84e-02 | ord  1.07  1.05
rw_t1000     kn1  hesdirk2 C f1     lev 32-64-128-256  L1 6.32e-02 3.35e-02 1.78e-02 | ord  0.92  0.91
rw_t1000     kn1  hesdirk2 C dens   lev 32-64-128-256  L1 7.95e-02 3.85e-02 1.88e-02 | ord  1.05  1.04
rw_t1000     kn1  hesdirk2 C velx   lev 32-64-128-256  L1 6.43e-02 3.36e-02 1.72e-02 | ord  0.94  0.96
rw_t1000     kn1  hesdirk2 C eint   lev 32-64-128-256  L1 7.92e-02 3.82e-02 1.85e-02 | ord  1.05  1.04
rw_t1000     knp  be       X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     knp  be       X f1     lev 32-64-128-256  L1 4.93e-03 1.47e-02 5.91e-04 | ord -1.58  4.64
rw_t1000     knp  be       X dens   lev 32-64-128-256  L1 8.73e-03 2.04e-03 5.31e-04 | ord  2.10  1.94
rw_t1000     knp  be       X velx   lev 32-64-128-256  L1 8.40e-03 2.03e-03 5.36e-04 | ord  2.05  1.92
rw_t1000     knp  be       X eint   lev 32-64-128-256  L1 8.69e-03 2.02e-03 5.30e-04 | ord  2.10  1.93
rw_t1000     knp  be       T e      lev 2-4-8-16  L1 3.30e-03 1.68e-03 8.50e-04 | ord  0.97  0.98
rw_t1000     knp  be       T f1     lev 2-4-8-16  L1 2.77e-02 1.69e-03 8.52e-04 | ord  4.03  0.99
rw_t1000     knp  be       T dens   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.41e-04 | ord  0.98  0.99
rw_t1000     knp  be       T velx   lev 2-4-8-16  L1 3.31e-03 1.68e-03 8.50e-04 | ord  0.97  0.99
rw_t1000     knp  be       T eint   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.43e-04 | ord  0.97  0.99
rw_t1000     knp  be       C e      lev 32-64-128-256  L1 2.49e-02 1.15e-02 5.58e-03 | ord  1.12  1.04
rw_t1000     knp  be       C f1     lev 32-64-128-256  L1 2.62e-02 1.31e-02 2.99e-02 | ord  1.00 -1.19
rw_t1000     knp  be       C dens   lev 32-64-128-256  L1 2.62e-02 1.17e-02 5.62e-03 | ord  1.16  1.06
rw_t1000     knp  be       C velx   lev 32-64-128-256  L1 2.51e-02 1.17e-02 5.67e-03 | ord  1.10  1.05
rw_t1000     knp  be       C eint   lev 32-64-128-256  L1 2.56e-02 1.16e-02 5.60e-03 | ord  1.14  1.05
rw_t1000     knp  hesdirk2 X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     knp  hesdirk2 X f1     lev 32-64-128-256  L1 5.75e-03 5.87e-03 5.91e-04 | ord -0.03  3.31
rw_t1000     knp  hesdirk2 X dens   lev 32-64-128-256  L1 8.76e-03 2.04e-03 5.27e-04 | ord  2.11  1.95
rw_t1000     knp  hesdirk2 X velx   lev 32-64-128-256  L1 8.44e-03 2.04e-03 5.37e-04 | ord  2.05  1.92
rw_t1000     knp  hesdirk2 X eint   lev 32-64-128-256  L1 8.72e-03 2.01e-03 5.27e-04 | ord  2.12  1.93
rw_t1000     knp  hesdirk2 T e      lev 2-4-8-16  L1 1.14e-02 5.66e-03 2.81e-03 | ord  1.01  1.01
rw_t1000     knp  hesdirk2 T f1     lev 2-4-8-16  L1 1.20e-02 5.54e-03 2.79e-03 | ord  1.12  0.99
rw_t1000     knp  hesdirk2 T dens   lev 2-4-8-16  L1 1.16e-02 5.78e-03 2.88e-03 | ord  1.01  1.00
rw_t1000     knp  hesdirk2 T velx   lev 2-4-8-16  L1 1.09e-02 5.54e-03 2.78e-03 | ord  0.98  0.99
rw_t1000     knp  hesdirk2 T eint   lev 2-4-8-16  L1 1.14e-02 5.68e-03 2.83e-03 | ord  1.01  1.00
rw_t1000     knp  hesdirk2 C e      lev 32-64-128-256  L1 7.99e-02 3.81e-02 1.84e-02 | ord  1.07  1.05
rw_t1000     knp  hesdirk2 C f1     lev 32-64-128-256  L1 6.32e-02 3.35e-02 1.78e-02 | ord  0.92  0.91
rw_t1000     knp  hesdirk2 C dens   lev 32-64-128-256  L1 7.95e-02 3.85e-02 1.88e-02 | ord  1.05  1.04
rw_t1000     knp  hesdirk2 C velx   lev 32-64-128-256  L1 6.43e-02 3.36e-02 1.72e-02 | ord  0.94  0.96
rw_t1000     knp  hesdirk2 C eint   lev 32-64-128-256  L1 7.92e-02 3.82e-02 1.85e-02 | ord  1.05  1.04
```

## RESULTS/atm_hopf.txt (raw)

```
cen n=  32  L1 1.398e-02  max|dE/E| 1.633e-01  max(tau<1) 1.633e-01  max(tau>1) 5.547e-02  top 1.484e-01  tau_cell top 4.0e-03 bot 2.8e+01
cen n=  64  L1 9.676e-03  max|dE/E| 1.675e-01  max(tau<1) 1.675e-01  max(tau>1) 5.224e-02  top 1.516e-01  tau_cell top 2.0e-03 bot 1.4e+01
cen n= 128  L1 8.664e-03  max|dE/E| 1.697e-01  max(tau<1) 1.697e-01  max(tau>1) 5.578e-02  top 1.531e-01  tau_cell top 1.0e-03 bot 7.0e+00
cen n= 256  L1 8.443e-03  max|dE/E| 1.708e-01  max(tau<1) 1.708e-01  max(tau>1) 5.544e-02  top 1.539e-01  tau_cell top 5.0e-04 bot 3.5e+00
cen n= 512  L1 8.396e-03  max|dE/E| 1.713e-01  max(tau<1) 1.713e-01  max(tau>1) 5.590e-02  top 1.543e-01  tau_cell top 2.5e-04 bot 1.7e+00
cen Hopf L1 orders 0.53 0.16 0.04 0.01
cen self m1_e 1.82e-03 4.59e-04 1.18e-04 3.86e-05 | ord 1.99 1.96 1.61
cen self m1_f1 9.17e-01 7.37e-01 1.12e+00 5.44e+00 | ord 0.32 -0.60 -2.28
hr n=  32  L1 5.472e-03  max|dE/E| 1.788e-02  max(tau<1) 1.788e-02  max(tau>1) 1.206e-02  top 8.650e-04  tau_cell top 4.0e-03 bot 2.8e+01
hr n=  64  L1 1.498e-03  max|dE/E| 1.722e-02  max(tau<1) 1.722e-02  max(tau>1) 1.325e-02  top 6.658e-04  tau_cell top 2.0e-03 bot 1.4e+01
hr n= 128  L1 1.946e-03  max|dE/E| 1.774e-02  max(tau<1) 1.774e-02  max(tau>1) 1.468e-02  top 8.608e-04  tau_cell top 1.0e-03 bot 7.0e+00
hr n= 256  L1 3.054e-03  max|dE/E| 1.956e-02  max(tau<1) 1.956e-02  max(tau>1) 1.303e-02  top 9.828e-04  tau_cell top 5.0e-04 bot 3.5e+00
hr n= 512  L1 4.236e-03  max|dE/E| 2.015e-02  max(tau<1) 2.015e-02  max(tau>1) 9.668e-03  top 8.396e-04  tau_cell top 2.5e-04 bot 1.7e+00
hr Hopf L1 orders 1.87 -0.38 -0.65 -0.47
hr self m1_e 2.11e-03 9.83e-04 8.98e-04 1.10e-03 | ord 1.10 0.13 -0.29
hr self m1_f1 1.27e+00 1.91e+01 5.60e+01 9.75e-01 | ord -3.91 -1.55 5.85
hrx0 n=  32  L1 5.936e-03  max|dE/E| 1.536e-02  max(tau<1) 1.536e-02  max(tau>1) 6.389e-03  top 8.197e-04  tau_cell top 4.0e-03 bot 2.8e+01
hrx0 n=  64  L1 1.684e-03  max|dE/E| 1.732e-02  max(tau<1) 1.732e-02  max(tau>1) 5.168e-03  top 2.821e-04  tau_cell top 2.0e-03 bot 1.4e+01
hrx0 n= 128  L1 6.658e-04  max|dE/E| 1.840e-02  max(tau<1) 1.840e-02  max(tau>1) 5.471e-03  top 8.283e-05  tau_cell top 1.0e-03 bot 7.0e+00
hrx0 n= 256  L1 4.304e-04  max|dE/E| 1.893e-02  max(tau<1) 1.893e-02  max(tau>1) 4.958e-03  top 1.520e-05  tau_cell top 5.0e-04 bot 3.5e+00
hrx0 n= 512  L1 3.714e-04  max|dE/E| 1.928e-02  max(tau<1) 1.928e-02  max(tau>1) 4.742e-03  top 3.278e-05  tau_cell top 2.5e-04 bot 1.7e+00
hrx0 Hopf L1 orders 1.82 1.34 0.63 0.21
hrx0 self m1_e 1.77e-03 4.41e-04 1.11e-04 3.91e-05 | ord 2.00 1.99 1.51
hrx0 self m1_f1 1.62e+00 1.02e+01 3.29e+01 9.05e-01 | ord -2.65 -1.69 5.18
sq n=  32  L1 6.655e-03  max|dE/E| 2.840e-02  max(tau<1) 2.840e-02  max(tau>1) 6.626e-03  top 1.458e-03  tau_cell top 4.0e-03 bot 2.8e+01
sq n=  64  L1 1.946e-03  max|dE/E| 2.379e-02  max(tau<1) 2.379e-02  max(tau>1) 1.969e-03  top 6.330e-04  tau_cell top 2.0e-03 bot 1.4e+01
sq n= 128  L1 7.760e-04  max|dE/E| 2.151e-02  max(tau<1) 2.151e-02  max(tau>1) 3.145e-03  top 2.557e-04  tau_cell top 1.0e-03 bot 7.0e+00
sq n= 256  L1 4.744e-04  max|dE/E| 2.051e-02  max(tau<1) 2.051e-02  max(tau>1) 3.686e-03  top 7.791e-05  tau_cell top 5.0e-04 bot 3.5e+00
sq n= 512  L1 3.874e-04  max|dE/E| 2.006e-02  max(tau<1) 2.006e-02  max(tau>1) 4.058e-03  top 4.134e-06  tau_cell top 2.5e-04 bot 1.7e+00
sq Hopf L1 orders 1.77 1.33 0.71 0.29
sq self m1_e 2.06e-03 5.52e-04 1.63e-04 6.14e-05 | ord 1.90 1.76 1.41
sq self m1_f1 1.47e+00 6.63e+00 1.56e+01 1.34e+00 | ord -2.17 -1.23 3.54
sqf n=  32  L1 6.655e-03  max|dE/E| 2.840e-02  max(tau<1) 2.840e-02  max(tau>1) 6.626e-03  top 1.458e-03  tau_cell top 4.0e-03 bot 2.8e+01
sqf n=  64  L1 1.946e-03  max|dE/E| 2.379e-02  max(tau<1) 2.379e-02  max(tau>1) 1.969e-03  top 6.330e-04  tau_cell top 2.0e-03 bot 1.4e+01
sqf n= 128  L1 7.760e-04  max|dE/E| 2.151e-02  max(tau<1) 2.151e-02  max(tau>1) 3.145e-03  top 2.557e-04  tau_cell top 1.0e-03 bot 7.0e+00
sqf n= 256  L1 4.744e-04  max|dE/E| 2.051e-02  max(tau<1) 2.051e-02  max(tau>1) 3.686e-03  top 7.791e-05  tau_cell top 5.0e-04 bot 3.5e+00
sqf n= 512  L1 3.874e-04  max|dE/E| 2.006e-02  max(tau<1) 2.006e-02  max(tau>1) 4.058e-03  top 4.134e-06  tau_cell top 2.5e-04 bot 1.7e+00
sqf Hopf L1 orders 1.77 1.33 0.71 0.29
sqf self m1_e 2.06e-03 5.52e-04 1.63e-04 6.14e-05 | ord 1.90 1.76 1.41
sqf self m1_f1 1.47e+00 6.63e+00 1.56e+01 1.34e+00 | ord -2.17 -1.23 3.54
sqv n=  32  L1 6.655e-03  max|dE/E| 2.840e-02  max(tau<1) 2.840e-02  max(tau>1) 6.626e-03  top 1.458e-03  tau_cell top 4.0e-03 bot 2.8e+01
sqv n=  64  L1 1.946e-03  max|dE/E| 2.379e-02  max(tau<1) 2.379e-02  max(tau>1) 1.969e-03  top 6.330e-04  tau_cell top 2.0e-03 bot 1.4e+01
sqv n= 128  L1 7.760e-04  max|dE/E| 2.151e-02  max(tau<1) 2.151e-02  max(tau>1) 3.145e-03  top 2.557e-04  tau_cell top 1.0e-03 bot 7.0e+00
sqv n= 256  L1 4.744e-04  max|dE/E| 2.051e-02  max(tau<1) 2.051e-02  max(tau>1) 3.686e-03  top 7.791e-05  tau_cell top 5.0e-04 bot 3.5e+00
sqv n= 512  L1 3.874e-04  max|dE/E| 2.006e-02  max(tau<1) 2.006e-02  max(tau>1) 4.058e-03  top 4.134e-06  tau_cell top 2.5e-04 bot 1.7e+00
sqv Hopf L1 orders 1.77 1.33 0.71 0.29
sqv self m1_e 2.06e-03 5.52e-04 1.63e-04 6.14e-05 | ord 1.90 1.76 1.41
sqv self m1_f1 1.47e+00 6.63e+00 1.56e+01 1.34e+00 | ord -2.17 -1.23 3.54
sqa n=  32  L1 6.655e-03  max|dE/E| 2.840e-02  max(tau<1) 2.840e-02  max(tau>1) 6.626e-03  top 1.458e-03  tau_cell top 4.0e-03 bot 2.8e+01
sqa n=  64  L1 1.946e-03  max|dE/E| 2.379e-02  max(tau<1) 2.379e-02  max(tau>1) 1.969e-03  top 6.330e-04  tau_cell top 2.0e-03 bot 1.4e+01
sqa n= 128  L1 7.760e-04  max|dE/E| 2.151e-02  max(tau<1) 2.151e-02  max(tau>1) 3.145e-03  top 2.557e-04  tau_cell top 1.0e-03 bot 7.0e+00
sqa n= 256  L1 4.744e-04  max|dE/E| 2.051e-02  max(tau<1) 2.051e-02  max(tau>1) 3.686e-03  top 7.791e-05  tau_cell top 5.0e-04 bot 3.5e+00
sqa n= 512  L1 3.874e-04  max|dE/E| 2.006e-02  max(tau<1) 2.006e-02  max(tau>1) 4.058e-03  top 4.134e-06  tau_cell top 2.5e-04 bot 1.7e+00
sqa Hopf L1 orders 1.77 1.33 0.71 0.29
sqa self m1_e 2.06e-03 5.52e-04 1.63e-04 6.14e-05 | ord 1.90 1.76 1.41
sqa self m1_f1 1.47e+00 6.63e+00 1.56e+01 1.34e+00 | ord -2.17 -1.23 3.54
aq n=  32  L1 6.712e-03  max|dE/E| 2.845e-02  max(tau<1) 2.845e-02  max(tau>1) 6.863e-03  top 1.434e-03  tau_cell top 4.0e-03 bot 2.8e+01
aq n=  64  L1 1.965e-03  max|dE/E| 2.356e-02  max(tau<1) 2.356e-02  max(tau>1) 2.115e-03  top 6.013e-04  tau_cell top 2.0e-03 bot 1.4e+01
aq n= 128  L1 7.818e-04  max|dE/E| 2.145e-02  max(tau<1) 2.145e-02  max(tau>1) 3.003e-03  top 2.343e-04  tau_cell top 1.0e-03 bot 7.0e+00
aq n= 256  L1 4.774e-04  max|dE/E| 2.050e-02  max(tau<1) 2.050e-02  max(tau>1) 3.645e-03  top 6.986e-05  tau_cell top 5.0e-04 bot 3.5e+00
aq n= 512  L1 3.888e-04  max|dE/E| 2.006e-02  max(tau<1) 2.006e-02  max(tau>1) 4.067e-03  top 4.807e-06  tau_cell top 2.5e-04 bot 1.7e+00
aq Hopf L1 orders 1.77 1.33 0.71 0.30
aq self m1_e 2.08e-03 5.67e-04 1.69e-04 6.41e-05 | ord 1.87 1.74 1.40
aq self m1_f1 1.51e+00 6.04e+00 1.49e+01 1.29e+00 | ord -2.00 -1.30 3.53
kn01 n=  32  L1 5.566e-03  max|dE/E| 1.788e-02  max(tau<1) 1.788e-02  max(tau>1) 1.206e-02  top 8.646e-04  tau_cell top 4.0e-03 bot 2.8e+01
kn01 n=  64  L1 1.620e-03  max|dE/E| 1.722e-02  max(tau<1) 1.722e-02  max(tau>1) 1.325e-02  top 6.650e-04  tau_cell top 2.0e-03 bot 1.4e+01
kn01 n= 128  L1 1.620e-03  max|dE/E| 1.774e-02  max(tau<1) 1.774e-02  max(tau>1) 1.468e-02  top 8.596e-04  tau_cell top 1.0e-03 bot 7.0e+00
kn01 n= 256  L1 2.380e-03  max|dE/E| 1.955e-02  max(tau<1) 1.955e-02  max(tau>1) 1.303e-02  top 9.814e-04  tau_cell top 5.0e-04 bot 3.5e+00
kn01 n= 512  L1 2.928e-03  max|dE/E| 2.015e-02  max(tau<1) 2.015e-02  max(tau>1) 9.676e-03  top 8.381e-04  tau_cell top 2.5e-04 bot 1.7e+00
kn01 Hopf L1 orders 1.78 -0.00 -0.56 -0.30
kn01 self m1_e 2.04e-03 8.45e-04 6.26e-04 5.77e-04 | ord 1.27 0.43 0.12
kn01 self m1_f1 1.27e+00 1.88e+01 5.45e+01 1.01e+00 | ord -3.89 -1.54 5.76
kn03 n=  32  L1 5.668e-03  max|dE/E| 1.790e-02  max(tau<1) 1.790e-02  max(tau>1) 1.204e-02  top 8.436e-04  tau_cell top 4.0e-03 bot 2.8e+01
kn03 n=  64  L1 1.655e-03  max|dE/E| 1.725e-02  max(tau<1) 1.725e-02  max(tau>1) 1.320e-02  top 6.295e-04  tau_cell top 2.0e-03 bot 1.4e+01
kn03 n= 128  L1 1.239e-03  max|dE/E| 1.766e-02  max(tau<1) 1.766e-02  max(tau>1) 1.467e-02  top 7.991e-04  tau_cell top 1.0e-03 bot 7.0e+00
kn03 n= 256  L1 1.602e-03  max|dE/E| 1.943e-02  max(tau<1) 1.943e-02  max(tau>1) 1.311e-02  top 8.950e-04  tau_cell top 5.0e-04 bot 3.5e+00
kn03 n= 512  L1 1.604e-03  max|dE/E| 2.002e-02  max(tau<1) 2.002e-02  max(tau>1) 9.965e-03  top 7.425e-04  tau_cell top 2.5e-04 bot 1.7e+00
kn03 Hopf L1 orders 1.78 0.42 -0.37 -0.00
kn03 self m1_e 1.95e-03 6.74e-04 3.22e-04 1.33e-04 | ord 1.53 1.07 1.27
kn03 self m1_f1 1.27e+00 1.73e+01 4.54e+01 1.27e+00 | ord -3.77 -1.39 5.16
kn1 n=  32  L1 5.780e-03  max|dE/E| 1.660e-02  max(tau<1) 1.660e-02  max(tau>1) 1.003e-02  top 6.573e-04  tau_cell top 4.0e-03 bot 2.8e+01
kn1 n=  64  L1 1.629e-03  max|dE/E| 1.489e-02  max(tau<1) 1.489e-02  max(tau>1) 9.763e-03  top 2.768e-04  tau_cell top 2.0e-03 bot 1.4e+01
kn1 n= 128  L1 7.466e-04  max|dE/E| 1.689e-02  max(tau<1) 1.689e-02  max(tau>1) 9.867e-03  top 2.050e-04  tau_cell top 1.0e-03 bot 7.0e+00
kn1 n= 256  L1 5.909e-04  max|dE/E| 1.832e-02  max(tau<1) 1.832e-02  max(tau>1) 8.089e-03  top 1.554e-04  tau_cell top 5.0e-04 bot 3.5e+00
kn1 n= 512  L1 5.226e-04  max|dE/E| 1.902e-02  max(tau<1) 1.902e-02  max(tau>1) 6.599e-03  top 7.354e-05  tau_cell top 2.5e-04 bot 1.7e+00
kn1 Hopf L1 orders 1.83 1.13 0.34 0.18
kn1 self m1_e 1.80e-03 4.70e-04 1.55e-04 1.05e-04 | ord 1.94 1.60 0.56
kn1 self m1_f1 1.46e+00 1.28e+01 1.45e+02 5.74e-01 | ord -3.13 -3.51 7.99
knp n=  32  L1 6.713e-03  max|dE/E| 2.845e-02  max(tau<1) 2.845e-02  max(tau>1) 6.859e-03  top 1.434e-03  tau_cell top 4.0e-03 bot 2.8e+01
knp n=  64  L1 1.964e-03  max|dE/E| 2.356e-02  max(tau<1) 2.356e-02  max(tau>1) 2.111e-03  top 6.017e-04  tau_cell top 2.0e-03 bot 1.4e+01
knp n= 128  L1 7.812e-04  max|dE/E| 2.145e-02  max(tau<1) 2.145e-02  max(tau>1) 3.007e-03  top 2.348e-04  tau_cell top 1.0e-03 bot 7.0e+00
knp n= 256  L1 4.768e-04  max|dE/E| 2.050e-02  max(tau<1) 2.050e-02  max(tau>1) 3.648e-03  top 7.035e-05  tau_cell top 5.0e-04 bot 3.5e+00
knp n= 512  L1 3.883e-04  max|dE/E| 2.006e-02  max(tau<1) 2.006e-02  max(tau>1) 4.069e-03  top 4.419e-06  tau_cell top 2.5e-04 bot 1.7e+00
knp Hopf L1 orders 1.77 1.33 0.71 0.30
knp self m1_e 2.08e-03 5.67e-04 1.69e-04 6.40e-05 | ord 1.87 1.74 1.40
knp self m1_f1 1.51e+00 6.01e+00 1.45e+01 1.33e+00 | ord -1.99 -1.27 3.44
```

## RESULTS/atm_self_regions.txt (raw)

```
cen
  tau 0-0.1      e 3.89e-03 1.87e-03 9.18e-04 4.56e-04 | ord  1.06  1.02  1.01
  tau 0.1-1      e 2.10e-03 1.09e-03 5.95e-04 3.03e-04 | ord  0.95  0.87  0.97
  tau 1-10       e 1.37e-03 2.62e-04 7.58e-05 4.31e-05 | ord  2.39  1.79  0.81
  tau 10-1e+09   e 2.28e-03 5.48e-04 1.25e-04 3.31e-05 | ord  2.06  2.13  1.91
  Kn  0-0.03     e 2.33e-03 5.65e-04 1.30e-04 3.81e-05 | ord  2.04  2.12  1.78
  Kn  0.03-0.3   e 1.99e-03 4.43e-04 9.17e-05 1.96e-05 | ord  2.17  2.27  2.23
  Kn  0.3-1e+09  e 1.80e-03 9.91e-04 5.22e-04 2.72e-04 | ord  0.86  0.92  0.94
hr
  tau 0-0.1      e 3.08e-03 2.57e-03 1.49e-03 4.60e-04 | ord  0.26  0.79  1.70
  tau 0.1-1      e 2.93e-03 2.69e-03 3.41e-03 2.83e-03 | ord  0.12 -0.34  0.27
  tau 1-10       e 2.39e-03 1.78e-03 1.70e-03 1.63e-03 | ord  0.42  0.06  0.06
  tau 10-1e+09   e 2.57e-03 1.07e-03 9.57e-04 1.28e-03 | ord  1.26  0.17 -0.42
  Kn  0-0.03     e 2.55e-03 9.63e-04 7.89e-04 1.09e-03 | ord  1.41  0.29 -0.46
  Kn  0.03-0.3   e 2.63e-03 1.51e-03 1.57e-03 1.78e-03 | ord  0.80 -0.06 -0.18
  Kn  0.3-1e+09  e 2.36e-03 2.14e-03 1.97e-03 1.98e-03 | ord  0.14  0.12 -0.01
hrx0
  tau 0-0.1      e 1.60e-03 9.28e-04 4.53e-04 2.83e-04 | ord  0.79  1.03  0.68
  tau 0.1-1      e 3.09e-03 1.56e-03 7.77e-04 4.57e-04 | ord  0.99  1.00  0.77
  tau 1-10       e 1.24e-03 2.59e-04 1.00e-04 8.01e-05 | ord  2.26  1.37  0.32
  tau 10-1e+09   e 2.25e-03 5.34e-04 1.19e-04 3.01e-05 | ord  2.07  2.17  1.98
  Kn  0-0.03     e 2.31e-03 5.57e-04 1.26e-04 3.59e-05 | ord  2.05  2.14  1.81
  Kn  0.03-0.3   e 1.89e-03 4.09e-04 8.38e-05 2.19e-05 | ord  2.21  2.29  1.94
  Kn  0.3-1e+09  e 1.69e-03 9.75e-04 5.20e-04 3.16e-04 | ord  0.80  0.91  0.72
sq
  tau 0-0.1      e 4.05e-03 1.92e-03 8.87e-04 4.09e-04 | ord  1.08  1.11  1.12
  tau 0.1-1      e 3.85e-03 1.52e-03 8.74e-04 4.88e-04 | ord  1.34  0.80  0.84
  tau 1-10       e 2.63e-03 7.84e-04 2.89e-04 1.30e-04 | ord  1.75  1.44  1.15
  tau 10-1e+09   e 2.43e-03 6.17e-04 1.61e-04 5.31e-05 | ord  1.98  1.94  1.60
  Kn  0-0.03     e 2.43e-03 6.10e-04 1.54e-04 5.12e-05 | ord  1.99  1.99  1.59
  Kn  0.03-0.3   e 2.46e-03 6.60e-04 1.99e-04 7.08e-05 | ord  1.90  1.73  1.49
  Kn  0.3-1e+09  e 3.50e-03 1.39e-03 7.06e-04 3.65e-04 | ord  1.34  0.97  0.95
sqf
  tau 0-0.1      e 4.05e-03 1.92e-03 8.87e-04 4.09e-04 | ord  1.08  1.11  1.12
  tau 0.1-1      e 3.85e-03 1.52e-03 8.74e-04 4.88e-04 | ord  1.34  0.80  0.84
  tau 1-10       e 2.63e-03 7.84e-04 2.89e-04 1.30e-04 | ord  1.75  1.44  1.15
  tau 10-1e+09   e 2.43e-03 6.17e-04 1.61e-04 5.31e-05 | ord  1.98  1.94  1.60
  Kn  0-0.03     e 2.43e-03 6.10e-04 1.54e-04 5.12e-05 | ord  1.99  1.99  1.59
  Kn  0.03-0.3   e 2.46e-03 6.60e-04 1.99e-04 7.08e-05 | ord  1.90  1.73  1.49
  Kn  0.3-1e+09  e 3.50e-03 1.39e-03 7.06e-04 3.65e-04 | ord  1.34  0.97  0.95
sqv
  tau 0-0.1      e 4.05e-03 1.92e-03 8.87e-04 4.09e-04 | ord  1.08  1.11  1.12
  tau 0.1-1      e 3.85e-03 1.52e-03 8.74e-04 4.88e-04 | ord  1.34  0.80  0.84
  tau 1-10       e 2.63e-03 7.84e-04 2.89e-04 1.30e-04 | ord  1.75  1.44  1.15
  tau 10-1e+09   e 2.43e-03 6.17e-04 1.61e-04 5.31e-05 | ord  1.98  1.94  1.60
  Kn  0-0.03     e 2.43e-03 6.10e-04 1.54e-04 5.12e-05 | ord  1.99  1.99  1.59
  Kn  0.03-0.3   e 2.46e-03 6.60e-04 1.99e-04 7.08e-05 | ord  1.90  1.73  1.49
  Kn  0.3-1e+09  e 3.50e-03 1.39e-03 7.06e-04 3.65e-04 | ord  1.34  0.97  0.95
sqa
  tau 0-0.1      e 4.05e-03 1.92e-03 8.87e-04 4.09e-04 | ord  1.08  1.11  1.12
  tau 0.1-1      e 3.85e-03 1.52e-03 8.74e-04 4.88e-04 | ord  1.34  0.80  0.84
  tau 1-10       e 2.63e-03 7.84e-04 2.89e-04 1.30e-04 | ord  1.75  1.44  1.15
  tau 10-1e+09   e 2.43e-03 6.17e-04 1.61e-04 5.31e-05 | ord  1.98  1.94  1.60
  Kn  0-0.03     e 2.43e-03 6.10e-04 1.54e-04 5.12e-05 | ord  1.99  1.99  1.59
  Kn  0.03-0.3   e 2.46e-03 6.60e-04 1.99e-04 7.08e-05 | ord  1.90  1.73  1.49
  Kn  0.3-1e+09  e 3.50e-03 1.39e-03 7.06e-04 3.65e-04 | ord  1.34  0.97  0.95
aq
  tau 0-0.1      e 4.08e-03 1.85e-03 8.49e-04 3.98e-04 | ord  1.14  1.12  1.09
  tau 0.1-1      e 4.32e-03 1.81e-03 9.22e-04 4.71e-04 | ord  1.26  0.97  0.97
  tau 1-10       e 2.71e-03 8.54e-04 3.23e-04 1.49e-04 | ord  1.67  1.40  1.11
  tau 10-1e+09   e 2.44e-03 6.25e-04 1.65e-04 5.51e-05 | ord  1.96  1.93  1.58
  Kn  0-0.03     e 2.44e-03 6.15e-04 1.56e-04 5.25e-05 | ord  1.99  1.98  1.58
  Kn  0.03-0.3   e 2.49e-03 6.84e-04 2.12e-04 7.78e-05 | ord  1.86  1.69  1.44
  Kn  0.3-1e+09  e 3.75e-03 1.53e-03 7.43e-04 3.71e-04 | ord  1.29  1.05  1.00
kn01
  tau 0-0.1      e 3.08e-03 2.57e-03 1.49e-03 4.59e-04 | ord  0.26  0.79  1.70
  tau 0.1-1      e 2.93e-03 2.69e-03 3.41e-03 2.83e-03 | ord  0.12 -0.34  0.27
  tau 1-10       e 2.36e-03 1.74e-03 1.62e-03 1.50e-03 | ord  0.44  0.10  0.11
  tau 10-1e+09   e 2.48e-03 8.86e-04 5.84e-04 5.60e-04 | ord  1.48  0.60  0.06
  Kn  0-0.03     e 2.46e-03 7.84e-04 4.35e-04 3.91e-04 | ord  1.65  0.85  0.15
  Kn  0.03-0.3   e 2.54e-03 1.33e-03 1.23e-03 1.15e-03 | ord  0.93  0.12  0.10
  Kn  0.3-1e+09  e 2.36e-03 2.14e-03 1.96e-03 1.97e-03 | ord  0.14  0.13 -0.01
kn03
  tau 0-0.1      e 3.06e-03 2.54e-03 1.46e-03 4.52e-04 | ord  0.27  0.80  1.69
  tau 0.1-1      e 2.92e-03 2.65e-03 3.32e-03 2.74e-03 | ord  0.14 -0.32  0.27
  tau 1-10       e 2.09e-03 1.23e-03 7.66e-04 5.09e-04 | ord  0.76  0.69  0.59
  tau 10-1e+09   e 2.39e-03 7.03e-04 2.59e-04 5.38e-05 | ord  1.76  1.44  2.27
  Kn  0-0.03     e 2.40e-03 6.65e-04 2.19e-04 5.16e-05 | ord  1.85  1.60  2.09
  Kn  0.03-0.3   e 2.32e-03 9.13e-04 4.87e-04 7.54e-05 | ord  1.35  0.91  2.69
  Kn  0.3-1e+09  e 2.32e-03 2.06e-03 1.86e-03 1.90e-03 | ord  0.17  0.14 -0.03
kn1
  tau 0-0.1      e 2.81e-03 2.11e-03 1.16e-03 5.42e-04 | ord  0.41  0.86  1.10
  tau 0.1-1      e 3.43e-03 2.94e-03 2.57e-03 1.75e-03 | ord  0.23  0.19  0.55
  tau 1-10       e 1.31e-03 2.91e-04 4.65e-04 5.25e-04 | ord  2.18 -0.68 -0.17
  tau 10-1e+09   e 2.27e-03 5.19e-04 8.36e-05 3.58e-05 | ord  2.13  2.64  1.22
  Kn  0-0.03     e 2.32e-03 5.47e-04 9.84e-05 1.79e-05 | ord  2.09  2.47  2.46
  Kn  0.03-0.3   e 1.95e-03 3.77e-04 1.08e-04 1.62e-04 | ord  2.37  1.80 -0.58
  Kn  0.3-1e+09  e 2.04e-03 1.88e-03 1.69e-03 1.22e-03 | ord  0.12  0.15  0.47
knp
  tau 0-0.1      e 4.08e-03 1.85e-03 8.49e-04 3.98e-04 | ord  1.14  1.12  1.09
  tau 0.1-1      e 4.32e-03 1.81e-03 9.21e-04 4.71e-04 | ord  1.26  0.97  0.97
  tau 1-10       e 2.72e-03 8.54e-04 3.23e-04 1.49e-04 | ord  1.67  1.40  1.12
  tau 10-1e+09   e 2.44e-03 6.25e-04 1.65e-04 5.51e-05 | ord  1.96  1.92  1.58
  Kn  0-0.03     e 2.44e-03 6.15e-04 1.57e-04 5.25e-05 | ord  1.99  1.97  1.58
  Kn  0.03-0.3   e 2.49e-03 6.84e-04 2.12e-04 7.75e-05 | ord  1.86  1.69  1.45
  Kn  0.3-1e+09  e 3.75e-03 1.53e-03 7.43e-04 3.70e-04 | ord  1.29  1.05  1.00
```

## RESULTS/arm_minus_cen.txt (raw)

```
pulse_k0.128  hr   hesdirk2 X |hr-cen| 1.43e-05 6.13e-05 2.56e-04 1.05e-03 4.23e-03
pulse_k0.128  hr   hesdirk2 C |hr-cen| 1.05e-04 2.29e-04 4.80e-04 9.82e-04 1.99e-03
pulse_k0.128  hr   hesdirk2 T |hr-cen| 7.83e-03 3.92e-03 1.96e-03 9.82e-04 4.91e-04
pulse_k0.128  hr   be       X |hr-cen| 4.92e-05 2.10e-04 8.77e-04 3.58e-03 1.44e-02
pulse_k0.128  hr   be       C |hr-cen| 3.73e-04 7.89e-04 1.65e-03 3.35e-03 6.76e-03
pulse_k0.128  hr   be       T |hr-cen| 2.67e-02 1.34e-02 6.70e-03 3.35e-03 1.68e-03
pulse_k0.128  sq   hesdirk2 X |sq-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sq   hesdirk2 C |sq-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sq   hesdirk2 T |sq-cen| 9.94e-15 1.28e-14 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sq   be       X |sq-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sq   be       C |sq-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sq   be       T |sq-cen| 6.61e-13 2.39e-14 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sqf  hesdirk2 X |sqf-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sqf  hesdirk2 C |sqf-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sqf  hesdirk2 T |sqf-cen| 9.94e-15 1.28e-14 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sqf  be       X |sqf-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sqf  be       C |sqf-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sqf  be       T |sqf-cen| 6.61e-13 2.39e-14 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sqv  hesdirk2 X |sqv-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sqv  hesdirk2 C |sqv-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sqv  hesdirk2 T |sqv-cen| 9.94e-15 1.28e-14 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sqv  be       X |sqv-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sqv  be       C |sqv-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sqv  be       T |sqv-cen| 6.61e-13 2.39e-14 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sqa  hesdirk2 X |sqa-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sqa  hesdirk2 C |sqa-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sqa  hesdirk2 T |sqa-cen| 9.94e-15 1.28e-14 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sqa  be       X |sqa-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sqa  be       C |sqa-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  sqa  be       T |sqa-cen| 6.61e-13 2.39e-14 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  aq   hesdirk2 X |aq-cen| 1.61e-05 6.73e-05 2.68e-04 1.07e-03 4.28e-03
pulse_k0.128  aq   hesdirk2 C |aq-cen| 1.21e-04 2.52e-04 5.03e-04 1.01e-03 2.01e-03
pulse_k0.128  aq   hesdirk2 T |aq-cen| 8.05e-03 4.02e-03 2.01e-03 1.01e-03 5.03e-04
pulse_k0.128  aq   be       X |aq-cen| 5.53e-05 2.30e-04 9.18e-04 3.66e-03 1.45e-02
pulse_k0.128  aq   be       C |aq-cen| 4.19e-04 8.63e-04 1.72e-03 3.43e-03 6.84e-03
pulse_k0.128  aq   be       T |aq-cen| 2.73e-02 1.37e-02 6.86e-03 3.43e-03 1.72e-03
pulse_k0.128  kn01 hesdirk2 X |kn01-cen| 1.47e-05 6.24e-05 2.60e-04 1.06e-03 4.29e-03
pulse_k0.128  kn01 hesdirk2 C |kn01-cen| 1.08e-04 2.33e-04 4.87e-04 9.96e-04 2.01e-03
pulse_k0.128  kn01 hesdirk2 T |kn01-cen| 7.95e-03 3.98e-03 1.99e-03 9.96e-04 4.98e-04
pulse_k0.128  kn01 be       X |kn01-cen| 5.05e-05 2.14e-04 8.89e-04 3.59e-03 1.44e-02
pulse_k0.128  kn01 be       C |kn01-cen| 3.82e-04 8.00e-04 1.66e-03 3.37e-03 6.78e-03
pulse_k0.128  kn01 be       T |kn01-cen| 2.69e-02 1.35e-02 6.72e-03 3.37e-03 1.69e-03
pulse_k0.128  kn03 hesdirk2 X |kn03-cen| 1.45e-05 6.16e-05 2.57e-04 1.05e-03 4.25e-03
pulse_k0.128  kn03 hesdirk2 C |kn03-cen| 1.07e-04 2.29e-04 4.81e-04 9.87e-04 2.00e-03
pulse_k0.128  kn03 hesdirk2 T |kn03-cen| 7.85e-03 3.93e-03 1.97e-03 9.87e-04 4.94e-04
pulse_k0.128  kn03 be       X |kn03-cen| 4.97e-05 2.11e-04 8.77e-04 3.56e-03 1.43e-02
pulse_k0.128  kn03 be       C |kn03-cen| 3.77e-04 7.86e-04 1.64e-03 3.34e-03 6.73e-03
pulse_k0.128  kn03 be       T |kn03-cen| 2.66e-02 1.33e-02 6.67e-03 3.34e-03 1.67e-03
pulse_k0.128  kn1  hesdirk2 X |kn1-cen| 1.41e-05 6.03e-05 2.52e-04 1.03e-03 4.17e-03
pulse_k0.128  kn1  hesdirk2 C |kn1-cen| 1.04e-04 2.24e-04 4.72e-04 9.66e-04 1.96e-03
pulse_k0.128  kn1  hesdirk2 T |kn1-cen| 7.69e-03 3.85e-03 1.93e-03 9.66e-04 4.83e-04
pulse_k0.128  kn1  be       X |kn1-cen| 4.85e-05 2.06e-04 8.62e-04 3.51e-03 1.41e-02
pulse_k0.128  kn1  be       C |kn1-cen| 3.68e-04 7.73e-04 1.62e-03 3.29e-03 6.64e-03
pulse_k0.128  kn1  be       T |kn1-cen| 2.60e-02 1.31e-02 6.57e-03 3.29e-03 1.65e-03
pulse_k0.128  knp  hesdirk2 X |knp-cen| 1.62e-05 6.74e-05 2.69e-04 1.08e-03 4.30e-03
pulse_k0.128  knp  hesdirk2 C |knp-cen| 1.21e-04 2.52e-04 5.04e-04 1.01e-03 2.02e-03
pulse_k0.128  knp  hesdirk2 T |knp-cen| 8.07e-03 4.03e-03 2.02e-03 1.01e-03 5.05e-04
pulse_k0.128  knp  be       X |knp-cen| 5.55e-05 2.30e-04 9.16e-04 3.65e-03 1.45e-02
pulse_k0.128  knp  be       C |knp-cen| 4.21e-04 8.57e-04 1.71e-03 3.42e-03 6.80e-03
pulse_k0.128  knp  be       T |knp-cen| 2.72e-02 1.36e-02 6.82e-03 3.42e-03 1.71e-03
pulse_k12.8   hr   hesdirk2 X |hr-cen| 3.38e-06 1.21e-05 4.96e-05 2.02e-04 8.20e-04
pulse_k12.8   hr   hesdirk2 C |hr-cen| 2.59e-05 4.55e-05 9.23e-05 1.90e-04 3.87e-04
pulse_k12.8   hr   hesdirk2 T |hr-cen| 1.42e-03 7.37e-04 3.76e-04 1.90e-04 9.55e-05
pulse_k12.8   hr   be       X |hr-cen| 1.12e-05 4.03e-05 1.66e-04 6.78e-04 2.74e-03
pulse_k12.8   hr   be       C |hr-cen| 7.07e-05 1.40e-04 3.02e-04 6.36e-04 1.31e-03
pulse_k12.8   hr   be       T |hr-cen| 4.17e-03 2.34e-03 1.24e-03 6.36e-04 3.23e-04
pulse_k12.8   sq   hesdirk2 X |sq-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sq   hesdirk2 C |sq-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sq   hesdirk2 T |sq-cen| 8.62e-15 7.14e-15 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sq   be       X |sq-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sq   be       C |sq-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sq   be       T |sq-cen| 7.92e-14 2.36e-14 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sqf  hesdirk2 X |sqf-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sqf  hesdirk2 C |sqf-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sqf  hesdirk2 T |sqf-cen| 8.62e-15 7.14e-15 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sqf  be       X |sqf-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sqf  be       C |sqf-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sqf  be       T |sqf-cen| 7.92e-14 2.36e-14 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sqv  hesdirk2 X |sqv-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sqv  hesdirk2 C |sqv-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sqv  hesdirk2 T |sqv-cen| 8.62e-15 7.14e-15 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sqv  be       X |sqv-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sqv  be       C |sqv-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sqv  be       T |sqv-cen| 7.92e-14 2.36e-14 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sqa  hesdirk2 X |sqa-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sqa  hesdirk2 C |sqa-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sqa  hesdirk2 T |sqa-cen| 8.62e-15 7.14e-15 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sqa  be       X |sqa-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sqa  be       C |sqa-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   sqa  be       T |sqa-cen| 7.92e-14 2.36e-14 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   aq   hesdirk2 X |aq-cen| 3.22e-06 1.35e-05 5.29e-05 2.11e-04 8.37e-04
pulse_k12.8   aq   hesdirk2 C |aq-cen| 2.33e-05 4.90e-05 9.82e-05 1.98e-04 3.95e-04
pulse_k12.8   aq   hesdirk2 T |aq-cen| 1.46e-03 7.62e-04 3.90e-04 1.98e-04 9.94e-05
pulse_k12.8   aq   be       X |aq-cen| 1.07e-05 4.53e-05 1.78e-04 7.06e-04 2.80e-03
pulse_k12.8   aq   be       C |aq-cen| 6.79e-05 1.58e-04 3.25e-04 6.64e-04 1.34e-03
pulse_k12.8   aq   be       T |aq-cen| 4.37e-03 2.44e-03 1.29e-03 6.64e-04 3.37e-04
pulse_k12.8   kn01 hesdirk2 X |kn01-cen| 2.31e-06 7.09e-06 2.76e-05 1.14e-04 4.74e-04
pulse_k12.8   kn01 hesdirk2 C |kn01-cen| 2.01e-05 2.69e-05 5.10e-05 1.07e-04 2.25e-04
pulse_k12.8   kn01 hesdirk2 T |kn01-cen| 7.67e-04 4.06e-04 2.10e-04 1.07e-04 5.42e-05
pulse_k12.8   kn01 be       X |kn01-cen| 7.53e-06 2.30e-05 9.14e-05 3.86e-04 1.62e-03
pulse_k12.8   kn01 be       C |kn01-cen| 5.13e-05 7.76e-05 1.65e-04 3.63e-04 7.68e-04
pulse_k12.8   kn01 be       T |kn01-cen| 2.43e-03 1.34e-03 7.06e-04 3.63e-04 1.84e-04
pulse_k12.8   kn03 hesdirk2 X |kn03-cen| 1.69e-06 4.47e-06 1.39e-05 5.90e-05 2.54e-04
pulse_k12.8   kn03 hesdirk2 C |kn03-cen| 1.64e-05 1.88e-05 2.63e-05 5.54e-05 1.22e-04
pulse_k12.8   kn03 hesdirk2 T |kn03-cen| 4.21e-04 2.11e-04 1.09e-04 5.54e-05 2.80e-05
pulse_k12.8   kn03 be       X |kn03-cen| 5.46e-06 1.42e-05 4.52e-05 1.99e-04 8.68e-04
pulse_k12.8   kn03 be       C |kn03-cen| 3.73e-05 4.99e-05 8.16e-05 1.88e-04 4.16e-04
pulse_k12.8   kn03 be       T |kn03-cen| 1.14e-03 6.67e-04 3.61e-04 1.88e-04 9.56e-05
pulse_k12.8   kn1  hesdirk2 X |kn1-cen| 4.20e-07 1.31e-06 4.15e-06 1.37e-05 4.90e-05
pulse_k12.8   kn1  hesdirk2 C |kn1-cen| 4.30e-06 6.06e-06 8.47e-06 1.28e-05 2.15e-05
pulse_k12.8   kn1  hesdirk2 T |kn1-cen| 1.71e-04 6.68e-05 2.82e-05 1.28e-05 6.04e-06
pulse_k12.8   kn1  be       X |kn1-cen| 1.30e-06 4.01e-06 1.22e-05 3.92e-05 1.38e-04
pulse_k12.8   kn1  be       C |kn1-cen| 8.35e-06 1.44e-05 2.22e-05 3.69e-05 6.69e-05
pulse_k12.8   kn1  be       T |kn1-cen| 2.53e-04 1.33e-04 7.08e-05 3.69e-05 1.89e-05
pulse_k12.8   knp  hesdirk2 X |knp-cen| 8.31e-07 4.65e-06 1.80e-05 7.05e-05 2.79e-04
pulse_k12.8   knp  hesdirk2 C |knp-cen| 6.35e-06 1.59e-05 3.28e-05 6.62e-05 1.34e-04
pulse_k12.8   knp  hesdirk2 T |knp-cen| 4.28e-04 2.39e-04 1.28e-04 6.62e-05 3.37e-05
pulse_k12.8   knp  be       X |knp-cen| 2.81e-06 1.58e-05 6.11e-05 2.40e-04 9.53e-04
pulse_k12.8   knp  be       C |knp-cen| 1.69e-05 5.31e-05 1.10e-04 2.26e-04 4.56e-04
pulse_k12.8   knp  be       T |knp-cen| 1.40e-03 8.08e-04 4.35e-04 2.26e-04 1.15e-04
pulse_k128    hr   hesdirk2 X |hr-cen|    -        -        -        -        -    
pulse_k128    hr   hesdirk2 C |hr-cen|    -        -        -        -        -    
pulse_k128    hr   hesdirk2 T |hr-cen|    -        -        -        -        -    
pulse_k128    hr   be       X |hr-cen|    -        -        -        -        -    
pulse_k128    hr   be       C |hr-cen|    -        -        -        -        -    
pulse_k128    hr   be       T |hr-cen|    -        -        -        -        -    
pulse_k128    sq   hesdirk2 X |sq-cen|    -        -        -        -        -    
pulse_k128    sq   hesdirk2 C |sq-cen|    -        -        -        -        -    
pulse_k128    sq   hesdirk2 T |sq-cen|    -        -        -        -        -    
pulse_k128    sq   be       X |sq-cen|    -        -        -        -        -    
pulse_k128    sq   be       C |sq-cen|    -        -        -        -        -    
pulse_k128    sq   be       T |sq-cen|    -        -        -        -        -    
pulse_k128    sqf  hesdirk2 X |sqf-cen|    -        -        -        -        -    
pulse_k128    sqf  hesdirk2 C |sqf-cen|    -        -        -        -        -    
pulse_k128    sqf  hesdirk2 T |sqf-cen|    -        -        -        -        -    
pulse_k128    sqf  be       X |sqf-cen|    -        -        -        -        -    
pulse_k128    sqf  be       C |sqf-cen|    -        -        -        -        -    
pulse_k128    sqf  be       T |sqf-cen|    -        -        -        -        -    
pulse_k128    sqv  hesdirk2 X |sqv-cen|    -        -        -        -        -    
pulse_k128    sqv  hesdirk2 C |sqv-cen|    -        -        -        -        -    
pulse_k128    sqv  hesdirk2 T |sqv-cen|    -        -        -        -        -    
pulse_k128    sqv  be       X |sqv-cen|    -        -        -        -        -    
pulse_k128    sqv  be       C |sqv-cen|    -        -        -        -        -    
pulse_k128    sqv  be       T |sqv-cen|    -        -        -        -        -    
pulse_k128    sqa  hesdirk2 X |sqa-cen|    -        -        -        -        -    
pulse_k128    sqa  hesdirk2 C |sqa-cen|    -        -        -        -        -    
pulse_k128    sqa  hesdirk2 T |sqa-cen|    -        -        -        -        -    
pulse_k128    sqa  be       X |sqa-cen|    -        -        -        -        -    
pulse_k128    sqa  be       C |sqa-cen|    -        -        -        -        -    
pulse_k128    sqa  be       T |sqa-cen|    -        -        -        -        -    
pulse_k128    aq   hesdirk2 X |aq-cen|    -        -        -        -        -    
pulse_k128    aq   hesdirk2 C |aq-cen|    -        -        -        -        -    
pulse_k128    aq   hesdirk2 T |aq-cen|    -        -        -        -        -    
pulse_k128    aq   be       X |aq-cen|    -        -        -        -        -    
pulse_k128    aq   be       C |aq-cen|    -        -        -        -        -    
pulse_k128    aq   be       T |aq-cen|    -        -        -        -        -    
pulse_k128    kn01 hesdirk2 X |kn01-cen|    -        -        -        -        -    
pulse_k128    kn01 hesdirk2 C |kn01-cen|    -        -        -        -        -    
pulse_k128    kn01 hesdirk2 T |kn01-cen|    -        -        -        -        -    
pulse_k128    kn01 be       X |kn01-cen|    -        -        -        -        -    
pulse_k128    kn01 be       C |kn01-cen|    -        -        -        -        -    
pulse_k128    kn01 be       T |kn01-cen|    -        -        -        -        -    
pulse_k128    kn03 hesdirk2 X |kn03-cen|    -        -        -        -        -    
pulse_k128    kn03 hesdirk2 C |kn03-cen|    -        -        -        -        -    
pulse_k128    kn03 hesdirk2 T |kn03-cen|    -        -        -        -        -    
pulse_k128    kn03 be       X |kn03-cen|    -        -        -        -        -    
pulse_k128    kn03 be       C |kn03-cen|    -        -        -        -        -    
pulse_k128    kn03 be       T |kn03-cen|    -        -        -        -        -    
pulse_k128    kn1  hesdirk2 X |kn1-cen|    -        -        -        -        -    
pulse_k128    kn1  hesdirk2 C |kn1-cen|    -        -        -        -        -    
pulse_k128    kn1  hesdirk2 T |kn1-cen|    -        -        -        -        -    
pulse_k128    kn1  be       X |kn1-cen|    -        -        -        -        -    
pulse_k128    kn1  be       C |kn1-cen|    -        -        -        -        -    
pulse_k128    kn1  be       T |kn1-cen|    -        -        -        -        -    
pulse_k128    knp  hesdirk2 X |knp-cen|    -        -        -        -        -    
pulse_k128    knp  hesdirk2 C |knp-cen|    -        -        -        -        -    
pulse_k128    knp  hesdirk2 T |knp-cen|    -        -        -        -        -    
pulse_k128    knp  be       X |knp-cen|    -        -        -        -        -    
pulse_k128    knp  be       C |knp-cen|    -        -        -        -        -    
pulse_k128    knp  be       T |knp-cen|    -        -        -        -        -    
pulse_k1280   hr   hesdirk2 X |hr-cen| 4.12e-06 8.16e-06 1.62e-05 3.24e-05 6.49e-05
pulse_k1280   hr   hesdirk2 C |hr-cen| 4.71e-06 9.11e-06 1.74e-05 3.21e-05 5.58e-05
pulse_k1280   hr   hesdirk2 T |hr-cen| 3.72e-05 3.63e-05 3.48e-05 3.21e-05 2.79e-05
pulse_k1280   hr   be       X |hr-cen| 4.58e-06 9.06e-06 1.80e-05 3.60e-05 7.21e-05
pulse_k1280   hr   be       C |hr-cen| 4.74e-06 9.35e-06 1.84e-05 3.59e-05 6.85e-05
pulse_k1280   hr   be       T |hr-cen| 3.74e-05 3.73e-05 3.68e-05 3.59e-05 3.42e-05
pulse_k1280   sq   hesdirk2 X |sq-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sq   hesdirk2 C |sq-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sq   hesdirk2 T |sq-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sq   be       X |sq-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sq   be       C |sq-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sq   be       T |sq-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sqf  hesdirk2 X |sqf-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sqf  hesdirk2 C |sqf-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sqf  hesdirk2 T |sqf-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sqf  be       X |sqf-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sqf  be       C |sqf-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sqf  be       T |sqf-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sqv  hesdirk2 X |sqv-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sqv  hesdirk2 C |sqv-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sqv  hesdirk2 T |sqv-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sqv  be       X |sqv-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sqv  be       C |sqv-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sqv  be       T |sqv-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sqa  hesdirk2 X |sqa-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sqa  hesdirk2 C |sqa-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sqa  hesdirk2 T |sqa-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sqa  be       X |sqa-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sqa  be       C |sqa-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   sqa  be       T |sqa-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   aq   hesdirk2 X |aq-cen| 1.16e-06 5.60e-07 2.62e-07 1.15e-07 9.76e-08
pulse_k1280   aq   hesdirk2 C |aq-cen| 1.33e-06 6.29e-07 2.86e-07 1.12e-07 4.22e-08
pulse_k1280   aq   hesdirk2 T |aq-cen| 3.33e-07 2.08e-07 1.46e-07 1.12e-07 8.87e-08
pulse_k1280   aq   be       X |aq-cen| 1.29e-06 6.21e-07 2.91e-07 1.27e-07 1.08e-07
pulse_k1280   aq   be       C |aq-cen| 1.34e-06 6.45e-07 3.02e-07 1.25e-07 5.16e-08
pulse_k1280   aq   be       T |aq-cen| 3.30e-07 2.12e-07 1.54e-07 1.25e-07 1.08e-07
pulse_k1280   kn01 hesdirk2 X |kn01-cen| 6.81e-11 1.96e-10 4.71e-10 1.03e-09 2.17e-09
pulse_k1280   kn01 hesdirk2 C |kn01-cen| 7.76e-11 2.19e-10 5.05e-10 1.02e-09 1.87e-09
pulse_k1280   kn01 hesdirk2 T |kn01-cen| 1.18e-09 1.16e-09 1.11e-09 1.02e-09 8.90e-10
pulse_k1280   kn01 be       X |kn01-cen| 7.53e-11 2.17e-10 5.20e-10 1.14e-09 2.39e-09
pulse_k1280   kn01 be       C |kn01-cen| 7.62e-11 2.21e-10 5.28e-10 1.14e-09 2.28e-09
pulse_k1280   kn01 be       T |kn01-cen| 1.15e-09 1.16e-09 1.16e-09 1.14e-09 1.09e-09
pulse_k1280   kn03 hesdirk2 X |kn03-cen| 8.50e-13 2.45e-12 5.83e-12 1.28e-11 2.68e-11
pulse_k1280   kn03 hesdirk2 C |kn03-cen| 9.60e-13 2.71e-12 6.23e-12 1.27e-11 2.31e-11
pulse_k1280   kn03 hesdirk2 T |kn03-cen| 1.46e-11 1.43e-11 1.37e-11 1.27e-11 1.10e-11
pulse_k1280   kn03 be       X |kn03-cen| 9.29e-13 2.67e-12 6.42e-12 1.41e-11 2.95e-11
pulse_k1280   kn03 be       C |kn03-cen| 9.42e-13 2.73e-12 6.53e-12 1.40e-11 2.81e-11
pulse_k1280   kn03 be       T |kn03-cen| 1.42e-11 1.44e-11 1.43e-11 1.40e-11 1.34e-11
pulse_k1280   kn1  hesdirk2 X |kn1-cen| 4.96e-14 5.63e-14 5.69e-14 1.08e-13 2.25e-13
pulse_k1280   kn1  hesdirk2 C |kn1-cen| 1.23e-14 2.28e-14 4.99e-14 1.06e-13 1.97e-13
pulse_k1280   kn1  hesdirk2 T |kn1-cen| 1.18e-13 1.16e-13 1.14e-13 1.06e-13 1.08e-13
pulse_k1280   kn1  be       X |kn1-cen| 1.91e-14 1.90e-14 5.20e-14 1.15e-13 2.37e-13
pulse_k1280   kn1  be       C |kn1-cen| 6.87e-15 2.20e-14 5.47e-14 1.15e-13 2.33e-13
pulse_k1280   kn1  be       T |kn1-cen| 1.16e-13 1.16e-13 1.17e-13 1.15e-13 1.14e-13
pulse_k1280   knp  hesdirk2 X |knp-cen| 2.35e-13 1.61e-13 9.21e-14 4.92e-14 4.95e-14
pulse_k1280   knp  hesdirk2 C |knp-cen| 2.51e-13 1.61e-13 8.74e-14 4.31e-14 3.35e-14
pulse_k1280   knp  hesdirk2 T |knp-cen| 1.30e-13 7.63e-14 5.72e-14 4.31e-14 4.64e-14
pulse_k1280   knp  be       X |knp-cen| 2.50e-13 1.52e-13 8.34e-14 3.96e-14 4.71e-14
pulse_k1280   knp  be       C |knp-cen| 2.47e-13 1.58e-13 9.05e-14 4.31e-14 1.99e-14
pulse_k1280   knp  be       T |knp-cen| 1.23e-13 7.73e-14 5.11e-14 4.31e-14 3.82e-14
pulse_k12800  hr   hesdirk2 X |hr-cen|    -        -        -        -        -    
pulse_k12800  hr   hesdirk2 C |hr-cen|    -        -        -        -        -    
pulse_k12800  hr   hesdirk2 T |hr-cen|    -        -        -        -        -    
pulse_k12800  hr   be       X |hr-cen|    -        -        -        -        -    
pulse_k12800  hr   be       C |hr-cen|    -        -        -        -        -    
pulse_k12800  hr   be       T |hr-cen|    -        -        -        -        -    
pulse_k12800  sq   hesdirk2 X |sq-cen|    -        -        -        -        -    
pulse_k12800  sq   hesdirk2 C |sq-cen|    -        -        -        -        -    
pulse_k12800  sq   hesdirk2 T |sq-cen|    -        -        -        -        -    
pulse_k12800  sq   be       X |sq-cen|    -        -        -        -        -    
pulse_k12800  sq   be       C |sq-cen|    -        -        -        -        -    
pulse_k12800  sq   be       T |sq-cen|    -        -        -        -        -    
pulse_k12800  sqf  hesdirk2 X |sqf-cen|    -        -        -        -        -    
pulse_k12800  sqf  hesdirk2 C |sqf-cen|    -        -        -        -        -    
pulse_k12800  sqf  hesdirk2 T |sqf-cen|    -        -        -        -        -    
pulse_k12800  sqf  be       X |sqf-cen|    -        -        -        -        -    
pulse_k12800  sqf  be       C |sqf-cen|    -        -        -        -        -    
pulse_k12800  sqf  be       T |sqf-cen|    -        -        -        -        -    
pulse_k12800  sqv  hesdirk2 X |sqv-cen|    -        -        -        -        -    
pulse_k12800  sqv  hesdirk2 C |sqv-cen|    -        -        -        -        -    
pulse_k12800  sqv  hesdirk2 T |sqv-cen|    -        -        -        -        -    
pulse_k12800  sqv  be       X |sqv-cen|    -        -        -        -        -    
pulse_k12800  sqv  be       C |sqv-cen|    -        -        -        -        -    
pulse_k12800  sqv  be       T |sqv-cen|    -        -        -        -        -    
pulse_k12800  sqa  hesdirk2 X |sqa-cen|    -        -        -        -        -    
pulse_k12800  sqa  hesdirk2 C |sqa-cen|    -        -        -        -        -    
pulse_k12800  sqa  hesdirk2 T |sqa-cen|    -        -        -        -        -    
pulse_k12800  sqa  be       X |sqa-cen|    -        -        -        -        -    
pulse_k12800  sqa  be       C |sqa-cen|    -        -        -        -        -    
pulse_k12800  sqa  be       T |sqa-cen|    -        -        -        -        -    
pulse_k12800  aq   hesdirk2 X |aq-cen|    -        -        -        -        -    
pulse_k12800  aq   hesdirk2 C |aq-cen|    -        -        -        -        -    
pulse_k12800  aq   hesdirk2 T |aq-cen|    -        -        -        -        -    
pulse_k12800  aq   be       X |aq-cen|    -        -        -        -        -    
pulse_k12800  aq   be       C |aq-cen|    -        -        -        -        -    
pulse_k12800  aq   be       T |aq-cen|    -        -        -        -        -    
pulse_k12800  kn01 hesdirk2 X |kn01-cen|    -        -        -        -        -    
pulse_k12800  kn01 hesdirk2 C |kn01-cen|    -        -        -        -        -    
pulse_k12800  kn01 hesdirk2 T |kn01-cen|    -        -        -        -        -    
pulse_k12800  kn01 be       X |kn01-cen|    -        -        -        -        -    
pulse_k12800  kn01 be       C |kn01-cen|    -        -        -        -        -    
pulse_k12800  kn01 be       T |kn01-cen|    -        -        -        -        -    
pulse_k12800  kn03 hesdirk2 X |kn03-cen|    -        -        -        -        -    
pulse_k12800  kn03 hesdirk2 C |kn03-cen|    -        -        -        -        -    
pulse_k12800  kn03 hesdirk2 T |kn03-cen|    -        -        -        -        -    
pulse_k12800  kn03 be       X |kn03-cen|    -        -        -        -        -    
pulse_k12800  kn03 be       C |kn03-cen|    -        -        -        -        -    
pulse_k12800  kn03 be       T |kn03-cen|    -        -        -        -        -    
pulse_k12800  kn1  hesdirk2 X |kn1-cen|    -        -        -        -        -    
pulse_k12800  kn1  hesdirk2 C |kn1-cen|    -        -        -        -        -    
pulse_k12800  kn1  hesdirk2 T |kn1-cen|    -        -        -        -        -    
pulse_k12800  kn1  be       X |kn1-cen|    -        -        -        -        -    
pulse_k12800  kn1  be       C |kn1-cen|    -        -        -        -        -    
pulse_k12800  kn1  be       T |kn1-cen|    -        -        -        -        -    
pulse_k12800  knp  hesdirk2 X |knp-cen|    -        -        -        -        -    
pulse_k12800  knp  hesdirk2 C |knp-cen|    -        -        -        -        -    
pulse_k12800  knp  hesdirk2 T |knp-cen|    -        -        -        -        -    
pulse_k12800  knp  be       X |knp-cen|    -        -        -        -        -    
pulse_k12800  knp  be       C |knp-cen|    -        -        -        -        -    
pulse_k12800  knp  be       T |knp-cen|    -        -        -        -        -    
rw_t10        hr   hesdirk2 X |hr-cen| 8.09e-06 1.05e-06 5.90e-05 3.03e-04 8.68e-04
rw_t10        hr   hesdirk2 C |hr-cen| 2.53e-05 1.36e-05 9.20e-05 3.83e-04 8.61e-04
rw_t10        hr   hesdirk2 T |hr-cen| 9.18e-05 7.76e-05 5.76e-05 3.73e-05 2.19e-05
rw_t10        hr   be       X |hr-cen| 8.36e-05 2.74e-04 9.04e-04 2.59e-03 5.17e-03
rw_t10        hr   be       C |hr-cen| 1.51e-03 2.68e-03 4.67e-03 6.04e-03 4.98e-03
rw_t10        hr   be       T |hr-cen| 4.58e-03 2.06e-03 8.48e-04 2.98e-04 7.85e-05
rw_t10        sq   hesdirk2 X |sq-cen| 7.35e-07 6.63e-07 4.30e-07 2.66e-07 6.03e-06
rw_t10        sq   hesdirk2 C |sq-cen| 3.19e-07 9.88e-08 9.45e-08 1.86e-07 5.73e-06
rw_t10        sq   hesdirk2 T |sq-cen| 3.90e-08 2.97e-08 3.75e-07 9.22e-07 1.30e-08
rw_t10        sq   be       X |sq-cen| 1.82e-07 2.96e-08 2.80e-08 3.12e-07 1.80e-04
rw_t10        sq   be       C |sq-cen| 2.02e-07 1.38e-04 1.69e-05 3.77e-06 1.64e-04
rw_t10        sq   be       T |sq-cen| 1.61e-05 4.59e-06 1.71e-08 3.75e-07 1.47e-07
rw_t10        sqf  hesdirk2 X |sqf-cen| 7.35e-07 6.63e-07 4.30e-07 2.66e-07 6.03e-06
rw_t10        sqf  hesdirk2 C |sqf-cen| 3.19e-07 9.88e-08 9.45e-08 1.86e-07 5.73e-06
rw_t10        sqf  hesdirk2 T |sqf-cen| 3.90e-08 2.97e-08 3.75e-07 9.22e-07 1.30e-08
rw_t10        sqf  be       X |sqf-cen| 1.82e-07 2.96e-08 2.80e-08 3.12e-07 1.80e-04
rw_t10        sqf  be       C |sqf-cen| 2.02e-07 1.38e-04 1.69e-05 3.77e-06 1.64e-04
rw_t10        sqf  be       T |sqf-cen| 1.61e-05 4.59e-06 1.71e-08 3.75e-07 1.47e-07
rw_t10        sqv  hesdirk2 X |sqv-cen| 7.35e-07 6.63e-07 4.30e-07 2.66e-07 6.03e-06
rw_t10        sqv  hesdirk2 C |sqv-cen| 3.19e-07 9.88e-08 9.45e-08 1.86e-07 5.73e-06
rw_t10        sqv  hesdirk2 T |sqv-cen| 3.90e-08 2.97e-08 3.75e-07 9.22e-07 1.30e-08
rw_t10        sqv  be       X |sqv-cen| 1.82e-07 2.96e-08 2.80e-08 3.12e-07 1.80e-04
rw_t10        sqv  be       C |sqv-cen| 2.02e-07 1.38e-04 1.69e-05 3.77e-06 1.64e-04
rw_t10        sqv  be       T |sqv-cen| 1.61e-05 4.59e-06 1.71e-08 3.75e-07 1.47e-07
rw_t10        sqa  hesdirk2 X |sqa-cen| 7.35e-07 6.63e-07 4.30e-07 2.66e-07 6.03e-06
rw_t10        sqa  hesdirk2 C |sqa-cen| 3.19e-07 9.88e-08 9.45e-08 1.86e-07 5.73e-06
rw_t10        sqa  hesdirk2 T |sqa-cen| 3.90e-08 2.97e-08 3.75e-07 9.22e-07 1.30e-08
rw_t10        sqa  be       X |sqa-cen| 1.82e-07 2.96e-08 2.80e-08 3.12e-07 1.80e-04
rw_t10        sqa  be       C |sqa-cen| 2.02e-07 1.38e-04 1.69e-05 3.77e-06 1.64e-04
rw_t10        sqa  be       T |sqa-cen| 1.61e-05 4.59e-06 1.71e-08 3.75e-07 1.47e-07
rw_t10        aq   hesdirk2 X |aq-cen| 9.98e-06 3.50e-05 1.27e-04 4.19e-04 1.01e-03
rw_t10        aq   hesdirk2 C |aq-cen| 1.67e-05 6.23e-05 2.05e-04 5.33e-04 1.01e-03
rw_t10        aq   hesdirk2 T |aq-cen| 2.05e-04 1.69e-04 1.24e-04 8.01e-05 4.69e-05
rw_t10        aq   be       X |aq-cen| 5.23e-05 2.12e-04 7.95e-04 2.44e-03 5.05e-03
rw_t10        aq   be       C |aq-cen| 1.44e-03 2.58e-03 4.55e-03 5.91e-03 4.86e-03
rw_t10        aq   be       T |aq-cen| 4.46e-03 1.94e-03 7.40e-04 2.12e-04 1.76e-05
rw_t10        kn01 hesdirk2 X |kn01-cen| 7.33e-07 6.13e-07 4.32e-07 3.22e-07 4.91e-06
rw_t10        kn01 hesdirk2 C |kn01-cen| 2.73e-07 8.34e-08 7.10e-08 9.16e-08 4.65e-06
rw_t10        kn01 hesdirk2 T |kn01-cen| 4.03e-08 2.70e-08 3.72e-07 8.19e-08 1.39e-08
rw_t10        kn01 be       X |kn01-cen| 3.23e-08 2.88e-08 3.19e-08 1.08e-06 6.51e-07
rw_t10        kn01 be       C |kn01-cen| 2.76e-08 1.39e-04 1.67e-05 9.26e-07 4.48e-07
rw_t10        kn01 be       T |kn01-cen| 1.63e-05 4.53e-06 2.90e-08 3.53e-07 1.34e-07
rw_t10        kn03 hesdirk2 X |kn03-cen| 7.33e-07 6.13e-07 4.32e-07 3.22e-07 4.91e-06
rw_t10        kn03 hesdirk2 C |kn03-cen| 2.73e-07 8.34e-08 7.10e-08 9.16e-08 4.65e-06
rw_t10        kn03 hesdirk2 T |kn03-cen| 4.03e-08 2.70e-08 3.72e-07 8.19e-08 1.39e-08
rw_t10        kn03 be       X |kn03-cen| 3.23e-08 2.88e-08 3.19e-08 1.08e-06 6.51e-07
rw_t10        kn03 be       C |kn03-cen| 2.76e-08 1.39e-04 1.67e-05 9.26e-07 4.48e-07
rw_t10        kn03 be       T |kn03-cen| 1.63e-05 4.53e-06 2.90e-08 3.53e-07 1.34e-07
rw_t10        kn1  hesdirk2 X |kn1-cen| 7.33e-07 6.13e-07 4.32e-07 3.22e-07 4.91e-06
rw_t10        kn1  hesdirk2 C |kn1-cen| 2.73e-07 8.34e-08 7.10e-08 9.16e-08 4.65e-06
rw_t10        kn1  hesdirk2 T |kn1-cen| 4.03e-08 2.70e-08 3.72e-07 8.19e-08 1.39e-08
rw_t10        kn1  be       X |kn1-cen| 3.23e-08 2.88e-08 3.19e-08 1.08e-06 6.51e-07
rw_t10        kn1  be       C |kn1-cen| 2.76e-08 1.39e-04 1.67e-05 9.26e-07 4.48e-07
rw_t10        kn1  be       T |kn1-cen| 1.63e-05 4.53e-06 2.90e-08 3.53e-07 1.34e-07
rw_t10        knp  hesdirk2 X |knp-cen| 7.33e-07 6.13e-07 4.32e-07 3.22e-07 4.91e-06
rw_t10        knp  hesdirk2 C |knp-cen| 2.73e-07 8.34e-08 7.10e-08 9.16e-08 4.65e-06
rw_t10        knp  hesdirk2 T |knp-cen| 4.03e-08 2.70e-08 3.72e-07 8.19e-08 1.39e-08
rw_t10        knp  be       X |knp-cen| 3.23e-08 2.88e-08 3.19e-08 1.08e-06 6.51e-07
rw_t10        knp  be       C |knp-cen| 2.76e-08 1.39e-04 1.67e-05 9.26e-07 4.48e-07
rw_t10        knp  be       T |knp-cen| 1.63e-05 4.53e-06 2.90e-08 3.53e-07 1.34e-07
rw_t1000      hr   hesdirk2 X |hr-cen| 8.19e-05 1.64e-04 3.28e-04 6.56e-04
rw_t1000      hr   hesdirk2 C |hr-cen| 8.09e-05 1.65e-04 3.32e-04 6.60e-04
rw_t1000      hr   hesdirk2 T |hr-cen| 3.30e-04 3.28e-04 3.22e-04 3.10e-04
rw_t1000      hr   be       X |hr-cen| 8.33e-05 1.67e-04 3.33e-04 6.66e-04
rw_t1000      hr   be       C |hr-cen| 8.41e-05 1.68e-04 3.34e-04 6.68e-04
rw_t1000      hr   be       T |hr-cen| 3.34e-04 3.33e-04 3.31e-04 3.27e-04
rw_t1000      sq   hesdirk2 X |sq-cen| 2.24e-07 2.24e-07 2.24e-07 2.24e-07
rw_t1000      sq   hesdirk2 C |sq-cen| 3.03e-10 1.90e-08 2.34e-09 3.60e-07
rw_t1000      sq   hesdirk2 T |sq-cen| 5.73e-07 2.18e-07 9.88e-08 4.58e-08
rw_t1000      sq   be       X |sq-cen| 5.35e-07 4.54e-07 1.53e-07 2.10e-08
rw_t1000      sq   be       C |sq-cen| 1.27e-09 4.97e-10 6.13e-09 1.75e-09
rw_t1000      sq   be       T |sq-cen| 2.85e-09 2.08e-07 2.07e-07 8.59e-08
rw_t1000      sqf  hesdirk2 X |sqf-cen| 2.24e-07 2.24e-07 2.24e-07 2.24e-07
rw_t1000      sqf  hesdirk2 C |sqf-cen| 3.03e-10 1.90e-08 2.34e-09 3.60e-07
rw_t1000      sqf  hesdirk2 T |sqf-cen| 5.73e-07 2.18e-07 9.88e-08 4.58e-08
rw_t1000      sqf  be       X |sqf-cen| 5.35e-07 4.54e-07 1.53e-07 2.10e-08
rw_t1000      sqf  be       C |sqf-cen| 1.27e-09 4.97e-10 6.13e-09 1.75e-09
rw_t1000      sqf  be       T |sqf-cen| 2.85e-09 2.08e-07 2.07e-07 8.59e-08
rw_t1000      sqv  hesdirk2 X |sqv-cen| 2.24e-07 2.24e-07 2.24e-07 2.24e-07
rw_t1000      sqv  hesdirk2 C |sqv-cen| 3.03e-10 1.90e-08 2.34e-09 3.60e-07
rw_t1000      sqv  hesdirk2 T |sqv-cen| 5.73e-07 2.18e-07 9.88e-08 4.58e-08
rw_t1000      sqv  be       X |sqv-cen| 5.35e-07 4.54e-07 1.53e-07 2.10e-08
rw_t1000      sqv  be       C |sqv-cen| 1.27e-09 4.97e-10 6.13e-09 1.75e-09
rw_t1000      sqv  be       T |sqv-cen| 2.85e-09 2.08e-07 2.07e-07 8.59e-08
rw_t1000      sqa  hesdirk2 X |sqa-cen| 2.24e-07 2.24e-07 2.24e-07 2.24e-07
rw_t1000      sqa  hesdirk2 C |sqa-cen| 3.03e-10 1.90e-08 2.34e-09 3.60e-07
rw_t1000      sqa  hesdirk2 T |sqa-cen| 5.73e-07 2.18e-07 9.88e-08 4.58e-08
rw_t1000      sqa  be       X |sqa-cen| 5.35e-07 4.54e-07 1.53e-07 2.10e-08
rw_t1000      sqa  be       C |sqa-cen| 1.27e-09 4.97e-10 6.13e-09 1.75e-09
rw_t1000      sqa  be       T |sqa-cen| 2.85e-09 2.08e-07 2.07e-07 8.59e-08
rw_t1000      aq   hesdirk2 X |aq-cen| 1.49e-06 6.99e-07 7.89e-07 2.84e-06
rw_t1000      aq   hesdirk2 C |aq-cen| 1.75e-06 1.33e-06 2.27e-06 4.44e-06
rw_t1000      aq   hesdirk2 T |aq-cen| 1.58e-06 7.78e-07 4.69e-07 3.82e-07
rw_t1000      aq   be       X |aq-cen| 1.75e-06 9.34e-07 8.76e-07 3.00e-06
rw_t1000      aq   be       C |aq-cen| 1.86e-06 1.48e-06 2.41e-06 4.67e-06
rw_t1000      aq   be       T |aq-cen| 1.56e-06 8.59e-07 3.74e-07 3.98e-07
rw_t1000      kn01 hesdirk2 X |kn01-cen| 2.24e-07 2.24e-07 2.24e-07 2.24e-07
rw_t1000      kn01 hesdirk2 C |kn01-cen| 3.03e-10 1.90e-08 2.34e-09 3.60e-07
rw_t1000      kn01 hesdirk2 T |kn01-cen| 5.73e-07 2.18e-07 9.88e-08 4.58e-08
rw_t1000      kn01 be       X |kn01-cen| 5.35e-07 4.54e-07 1.53e-07 2.10e-08
rw_t1000      kn01 be       C |kn01-cen| 1.27e-09 4.97e-10 6.13e-09 1.75e-09
rw_t1000      kn01 be       T |kn01-cen| 2.85e-09 2.08e-07 2.07e-07 8.59e-08
rw_t1000      kn03 hesdirk2 X |kn03-cen| 2.24e-07 2.24e-07 2.24e-07 2.24e-07
rw_t1000      kn03 hesdirk2 C |kn03-cen| 3.03e-10 1.90e-08 2.34e-09 3.60e-07
rw_t1000      kn03 hesdirk2 T |kn03-cen| 5.73e-07 2.18e-07 9.88e-08 4.58e-08
rw_t1000      kn03 be       X |kn03-cen| 5.35e-07 4.54e-07 1.53e-07 2.10e-08
rw_t1000      kn03 be       C |kn03-cen| 1.27e-09 4.97e-10 6.13e-09 1.75e-09
rw_t1000      kn03 be       T |kn03-cen| 2.85e-09 2.08e-07 2.07e-07 8.59e-08
rw_t1000      kn1  hesdirk2 X |kn1-cen| 2.24e-07 2.24e-07 2.24e-07 2.24e-07
rw_t1000      kn1  hesdirk2 C |kn1-cen| 3.03e-10 1.90e-08 2.34e-09 3.60e-07
rw_t1000      kn1  hesdirk2 T |kn1-cen| 5.73e-07 2.18e-07 9.88e-08 4.58e-08
rw_t1000      kn1  be       X |kn1-cen| 5.35e-07 4.54e-07 1.53e-07 2.10e-08
rw_t1000      kn1  be       C |kn1-cen| 1.27e-09 4.97e-10 6.13e-09 1.75e-09
rw_t1000      kn1  be       T |kn1-cen| 2.85e-09 2.08e-07 2.07e-07 8.59e-08
rw_t1000      knp  hesdirk2 X |knp-cen| 2.24e-07 2.24e-07 2.24e-07 2.24e-07
rw_t1000      knp  hesdirk2 C |knp-cen| 3.03e-10 1.90e-08 2.34e-09 3.60e-07
rw_t1000      knp  hesdirk2 T |knp-cen| 5.73e-07 2.18e-07 9.88e-08 4.58e-08
rw_t1000      knp  be       X |knp-cen| 5.35e-07 4.54e-07 1.53e-07 2.10e-08
rw_t1000      knp  be       C |knp-cen| 1.27e-09 4.97e-10 6.13e-09 1.75e-09
rw_t1000      knp  be       T |knp-cen| 2.85e-09 2.08e-07 2.07e-07 8.59e-08
```

## RESULTS/gates.txt (raw)

```
g1|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/g1_aq_1e4|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/fix_atm2d.athinput|rad_m1/closure=vet_sc rad_m1/vet_tensor=full rad_m1/implicit_flux_faces=all rad_m1/implicit_flux_beam=halfrange rad_m1/implicit_flux=blend rad_m1/implicit_blend=idort_f rad_m1/implicit_blend_tau0=1.0 rad_m1/implicit_blend_alpha=1.0 rad_m1/implicit_blend_flo=0.3 rad_m1/implicit_blend_fhi=0.6 rad_m1/implicit_blend_xthin=30 rad_m1/implicit_hr_recon=plm rad_m1/implicit_hr_recon_qs=0.1 rad_m1/implicit_cfl=1e4 time/tlim=2000 output1/dt=2000 output2/dt=2000: line 3: 1807235 Segmentation fault      (core dumped) OMP_NUM_THREADS=1 nice -n 10 $ATHENA_CPU -i $IN -d $R rad_m1/transport=implicit time/cfl_number=1e6 $K > run.log 2>&1
g3|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/g3_sq_c1|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/fix_marsh2d.athinput|rad_m1/closure=vet_sc rad_m1/vet_tensor=full rad_m1/implicit_flux_faces=all rad_m1/implicit_flux_beam=halfrange rad_m1/implicit_flux=blend rad_m1/implicit_blend=idort_f rad_m1/implicit_blend_tau0=1.0 rad_m1/implicit_blend_alpha=1.0 rad_m1/implicit_blend_flo=0.3 rad_m1/implicit_blend_fhi=0.6 rad_m1/implicit_blend_xthin_mode=steep rad_m1/implicit_blend_xthin=60 rad_m1/implicit_hr_recon=plm rad_m1/implicit_hr_recon_qs=0.1 rad_m1/implicit_cfl=1 rad_m1/implicit_bc_x1min=marshak rad_m1/implicit_ebath_x1min=1.0e-10 rad_m1/implicit_bc_x1max=marshak: line 3: 1807225 Segmentation fault      (core dumped) OMP_NUM_THREADS=1 nice -n 10 $ATHENA_CPU -i $IN -d $R rad_m1/transport=implicit time/cfl_number=1e6 $K > run.log 2>&1
g1|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/g1_knp_1e4|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/fix_atm2d.athinput|rad_m1/closure=vet_sc rad_m1/vet_tensor=full rad_m1/implicit_flux_faces=all rad_m1/implicit_flux_beam=halfrange rad_m1/implicit_flux=blend rad_m1/implicit_blend=idort_f rad_m1/implicit_blend_tau0=1.0 rad_m1/implicit_blend_alpha=1.0 rad_m1/implicit_blend_flo=0.3 rad_m1/implicit_blend_fhi=0.6 rad_m1/implicit_blend_xthin=30 rad_m1/implicit_blend_xthin_mode=kn rad_m1/implicit_blend_xthin_r0=0.3 rad_m1/implicit_hr_recon=plm rad_m1/implicit_hr_recon_qs=0.1 rad_m1/implicit_cfl=1e4 time/tlim=2000 output1/dt=2000 output2/dt=2000: line 3: 1807281 Segmentation fault      (core dumped) OMP_NUM_THREADS=1 nice -n 10 $ATHENA_CPU -i $IN -d $R rad_m1/transport=implicit time/cfl_number=1e6 $K > run.log 2>&1
g3|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/g3_sq_c10|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/fix_marsh2d.athinput|rad_m1/closure=vet_sc rad_m1/vet_tensor=full rad_m1/implicit_flux_faces=all rad_m1/implicit_flux_beam=halfrange rad_m1/implicit_flux=blend rad_m1/implicit_blend=idort_f rad_m1/implicit_blend_tau0=1.0 rad_m1/implicit_blend_alpha=1.0 rad_m1/implicit_blend_flo=0.3 rad_m1/implicit_blend_fhi=0.6 rad_m1/implicit_blend_xthin_mode=steep rad_m1/implicit_blend_xthin=60 rad_m1/implicit_hr_recon=plm rad_m1/implicit_hr_recon_qs=0.1 rad_m1/implicit_cfl=10 rad_m1/implicit_bc_x1min=marshak rad_m1/implicit_ebath_x1min=1.0e-10 rad_m1/implicit_bc_x1max=marshak: line 3: 1807226 Segmentation fault      (core dumped) OMP_NUM_THREADS=1 nice -n 10 $ATHENA_CPU -i $IN -d $R rad_m1/transport=implicit time/cfl_number=1e6 $K > run.log 2>&1
g1|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/g1_sq_1e4|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/fix_atm2d.athinput|rad_m1/closure=vet_sc rad_m1/vet_tensor=full rad_m1/implicit_flux_faces=all rad_m1/implicit_flux_beam=halfrange rad_m1/implicit_flux=blend rad_m1/implicit_blend=idort_f rad_m1/implicit_blend_tau0=1.0 rad_m1/implicit_blend_alpha=1.0 rad_m1/implicit_blend_flo=0.3 rad_m1/implicit_blend_fhi=0.6 rad_m1/implicit_blend_xthin_mode=steep rad_m1/implicit_blend_xthin=60 rad_m1/implicit_hr_recon=plm rad_m1/implicit_hr_recon_qs=0.1 rad_m1/implicit_cfl=1e4 time/tlim=2000 output1/dt=2000 output2/dt=2000: line 3: 1807209 Segmentation fault      (core dumped) OMP_NUM_THREADS=1 nice -n 10 $ATHENA_CPU -i $IN -d $R rad_m1/transport=implicit time/cfl_number=1e6 $K > run.log 2>&1
g3|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/g3_knp_c1|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/fix_marsh2d.athinput|rad_m1/closure=vet_sc rad_m1/vet_tensor=full rad_m1/implicit_flux_faces=all rad_m1/implicit_flux_beam=halfrange rad_m1/implicit_flux=blend rad_m1/implicit_blend=idort_f rad_m1/implicit_blend_tau0=1.0 rad_m1/implicit_blend_alpha=1.0 rad_m1/implicit_blend_flo=0.3 rad_m1/implicit_blend_fhi=0.6 rad_m1/implicit_blend_xthin=30 rad_m1/implicit_blend_xthin_mode=kn rad_m1/implicit_blend_xthin_r0=0.3 rad_m1/implicit_hr_recon=plm rad_m1/implicit_hr_recon_qs=0.1 rad_m1/implicit_cfl=1 rad_m1/implicit_bc_x1min=marshak rad_m1/implicit_ebath_x1min=1.0e-10 rad_m1/implicit_bc_x1max=marshak: line 3: 1807289 Segmentation fault      (core dumped) OMP_NUM_THREADS=1 nice -n 10 $ATHENA_CPU -i $IN -d $R rad_m1/transport=implicit time/cfl_number=1e6 $K > run.log 2>&1
g1|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/g1_aq_1e2|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/fix_atm2d.athinput|rad_m1/closure=vet_sc rad_m1/vet_tensor=full rad_m1/implicit_flux_faces=all rad_m1/implicit_flux_beam=halfrange rad_m1/implicit_flux=blend rad_m1/implicit_blend=idort_f rad_m1/implicit_blend_tau0=1.0 rad_m1/implicit_blend_alpha=1.0 rad_m1/implicit_blend_flo=0.3 rad_m1/implicit_blend_fhi=0.6 rad_m1/implicit_blend_xthin=30 rad_m1/implicit_hr_recon=plm rad_m1/implicit_hr_recon_qs=0.1 rad_m1/implicit_cfl=1e2 time/tlim=2000 output1/dt=2000 output2/dt=2000: line 3: 1807232 Segmentation fault      (core dumped) OMP_NUM_THREADS=1 nice -n 10 $ATHENA_CPU -i $IN -d $R rad_m1/transport=implicit time/cfl_number=1e6 $K > run.log 2>&1
g1|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/g1_sq_1e2|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/fix_atm2d.athinput|rad_m1/closure=vet_sc rad_m1/vet_tensor=full rad_m1/implicit_flux_faces=all rad_m1/implicit_flux_beam=halfrange rad_m1/implicit_flux=blend rad_m1/implicit_blend=idort_f rad_m1/implicit_blend_tau0=1.0 rad_m1/implicit_blend_alpha=1.0 rad_m1/implicit_blend_flo=0.3 rad_m1/implicit_blend_fhi=0.6 rad_m1/implicit_blend_xthin_mode=steep rad_m1/implicit_blend_xthin=60 rad_m1/implicit_hr_recon=plm rad_m1/implicit_hr_recon_qs=0.1 rad_m1/implicit_cfl=1e2 time/tlim=2000 output1/dt=2000 output2/dt=2000: line 3: 1807206 Segmentation fault      (core dumped) OMP_NUM_THREADS=1 nice -n 10 $ATHENA_CPU -i $IN -d $R rad_m1/transport=implicit time/cfl_number=1e6 $K > run.log 2>&1
g1|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/g1_knp_1e2|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/fix_atm2d.athinput|rad_m1/closure=vet_sc rad_m1/vet_tensor=full rad_m1/implicit_flux_faces=all rad_m1/implicit_flux_beam=halfrange rad_m1/implicit_flux=blend rad_m1/implicit_blend=idort_f rad_m1/implicit_blend_tau0=1.0 rad_m1/implicit_blend_alpha=1.0 rad_m1/implicit_blend_flo=0.3 rad_m1/implicit_blend_fhi=0.6 rad_m1/implicit_blend_xthin=30 rad_m1/implicit_blend_xthin_mode=kn rad_m1/implicit_blend_xthin_r0=0.3 rad_m1/implicit_hr_recon=plm rad_m1/implicit_hr_recon_qs=0.1 rad_m1/implicit_cfl=1e2 time/tlim=2000 output1/dt=2000 output2/dt=2000: line 3: 1807283 Segmentation fault      (core dumped) OMP_NUM_THREADS=1 nice -n 10 $ATHENA_CPU -i $IN -d $R rad_m1/transport=implicit time/cfl_number=1e6 $K > run.log 2>&1
g3|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/g3_knp_c10|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/fix_marsh2d.athinput|rad_m1/closure=vet_sc rad_m1/vet_tensor=full rad_m1/implicit_flux_faces=all rad_m1/implicit_flux_beam=halfrange rad_m1/implicit_flux=blend rad_m1/implicit_blend=idort_f rad_m1/implicit_blend_tau0=1.0 rad_m1/implicit_blend_alpha=1.0 rad_m1/implicit_blend_flo=0.3 rad_m1/implicit_blend_fhi=0.6 rad_m1/implicit_blend_xthin=30 rad_m1/implicit_blend_xthin_mode=kn rad_m1/implicit_blend_xthin_r0=0.3 rad_m1/implicit_hr_recon=plm rad_m1/implicit_hr_recon_qs=0.1 rad_m1/implicit_cfl=10 rad_m1/implicit_bc_x1min=marshak rad_m1/implicit_ebath_x1min=1.0e-10 rad_m1/implicit_bc_x1max=marshak: line 3: 1807290 Segmentation fault      (core dumped) OMP_NUM_THREADS=1 nice -n 10 $ATHENA_CPU -i $IN -d $R rad_m1/transport=implicit time/cfl_number=1e6 $K > run.log 2>&1
g3|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/g3_aq_c1|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/fix_marsh2d.athinput|rad_m1/closure=vet_sc rad_m1/vet_tensor=full rad_m1/implicit_flux_faces=all rad_m1/implicit_flux_beam=halfrange rad_m1/implicit_flux=blend rad_m1/implicit_blend=idort_f rad_m1/implicit_blend_tau0=1.0 rad_m1/implicit_blend_alpha=1.0 rad_m1/implicit_blend_flo=0.3 rad_m1/implicit_blend_fhi=0.6 rad_m1/implicit_blend_xthin=30 rad_m1/implicit_hr_recon=plm rad_m1/implicit_hr_recon_qs=0.1 rad_m1/implicit_cfl=1 rad_m1/implicit_bc_x1min=marshak rad_m1/implicit_ebath_x1min=1.0e-10 rad_m1/implicit_bc_x1max=marshak: line 3: 1807254 Segmentation fault      (core dumped) OMP_NUM_THREADS=1 nice -n 10 $ATHENA_CPU -i $IN -d $R rad_m1/transport=implicit time/cfl_number=1e6 $K > run.log 2>&1
g3|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/g3_aq_c10|/resnick/groups/carnegie_poc/jingze/arms_1010/run/gates/fix_marsh2d.athinput|rad_m1/closure=vet_sc rad_m1/vet_tensor=full rad_m1/implicit_flux_faces=all rad_m1/implicit_flux_beam=halfrange rad_m1/implicit_flux=blend rad_m1/implicit_blend=idort_f rad_m1/implicit_blend_tau0=1.0 rad_m1/implicit_blend_alpha=1.0 rad_m1/implicit_blend_flo=0.3 rad_m1/implicit_blend_fhi=0.6 rad_m1/implicit_blend_xthin=30 rad_m1/implicit_hr_recon=plm rad_m1/implicit_hr_recon_qs=0.1 rad_m1/implicit_cfl=10 rad_m1/implicit_bc_x1min=marshak rad_m1/implicit_ebath_x1min=1.0e-10 rad_m1/implicit_bc_x1max=marshak: line 3: 1807251 Segmentation fault      (core dumped) OMP_NUM_THREADS=1 nice -n 10 $ATHENA_CPU -i $IN -d $R rad_m1/transport=implicit time/cfl_number=1e6 $K > run.log 2>&1
g1_knp_1e4 rc=139  | HOPF: max|dE/E| 1.540e-01 top5 1.540e-01 L1 3.847e-03 (E_top 2.0015 exact 1.7344, tau_top-cell 5.09e-04) 
g1_aq_1e4 rc=139  | HOPF: max|dE/E| 1.540e-01 top5 1.540e-01 L1 3.847e-03 (E_top 2.0015 exact 1.7344, tau_top-cell 5.09e-04) 
g1_aq_1e2 rc=139  | HOPF: max|dE/E| 1.540e-01 top5 1.540e-01 L1 3.847e-03 (E_top 2.0015 exact 1.7344, tau_top-cell 5.09e-04) 
g1_knp_1e2 rc=139  | HOPF: max|dE/E| 1.540e-01 top5 1.540e-01 L1 3.847e-03 (E_top 2.0015 exact 1.7344, tau_top-cell 5.09e-04) 
g1_sq_1e2 rc=139  | HOPF: max|dE/E| 1.540e-01 top5 1.540e-01 L1 3.847e-03 (E_top 2.0015 exact 1.7344, tau_top-cell 5.09e-04) 
g1_sq_1e4 rc=139  | HOPF: max|dE/E| 1.540e-01 top5 1.540e-01 L1 3.847e-03 (E_top 2.0015 exact 1.7344, tau_top-cell 5.09e-04) 
g3_sq_c1 rc=139  | FAIL T6 marshak: L1(E)=1.0000 L1(material)=1.0000 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/arms_1010/repo/tests_m1/runs_3a/t6_ref_sn.txt 
g3_knp_c1 rc=139  | FAIL T6 marshak: L1(E)=1.0000 L1(material)=1.0000 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/arms_1010/repo/tests_m1/runs_3a/t6_ref_sn.txt 
g3_aq_c10 rc=139  | FAIL T6 marshak: L1(E)=1.0000 L1(material)=1.0000 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/arms_1010/repo/tests_m1/runs_3a/t6_ref_sn.txt 
g3_knp_c10 rc=139  | FAIL T6 marshak: L1(E)=1.0000 L1(material)=1.0000 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/arms_1010/repo/tests_m1/runs_3a/t6_ref_sn.txt 
g3_aq_c1 rc=139  | FAIL T6 marshak: L1(E)=1.0000 L1(material)=1.0000 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/arms_1010/repo/tests_m1/runs_3a/t6_ref_sn.txt 
g3_sq_c10 rc=139  | FAIL T6 marshak: L1(E)=1.0000 L1(material)=1.0000 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/arms_1010/repo/tests_m1/runs_3a/t6_ref_sn.txt 
g1_kn03_1e4 rc=0 Picard iterations mean=2.154902e+01 | HOPF: max|dE/E| 1.766e-02 top5 5.344e-03 L1 1.239e-03 (E_top 1.7358 exact 1.7344, tau_top-cell 5.09e-04) 
g1_hr_1e4 rc=0 Picard iterations mean=2.150980e+01 | HOPF: max|dE/E| 1.774e-02 top5 5.409e-03 L1 1.946e-03 (E_top 1.7359 exact 1.7344, tau_top-cell 5.09e-04) 
g5_hr_k12800 rc=0 Picard iterations mean=1.023529e+00 | PASS T3 pulse: d(sigma^2)/dt=5.208163e-05 2D=5.208333e-05 ratio=0.999967 (target 1.000000 +- 0.02) 
g5_hr_k1280 rc=0 Picard iterations mean=1.827451e+00 | PASS T3 pulse: d(sigma^2)/dt=5.295020e-04 2D=5.208333e-04 ratio=1.016644 (target 1.000000 +- 0.02) 
g5_kn03_k128 rc=0 Picard iterations mean=2.000000e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.780119e-03 2D=5.208333e-03 ratio=0.533783 (target 1.000000 +- 0.02) 
g5_kn03_k1280 rc=0 Picard iterations mean=1.827451e+00 | PASS T3 pulse: d(sigma^2)/dt=5.295065e-04 2D=5.208333e-04 ratio=1.016652 (target 1.000000 +- 0.02) 
g5_sq_k12800 rc=0 Picard iterations mean=1.023529e+00 | PASS T3 pulse: d(sigma^2)/dt=5.208134e-05 2D=5.208333e-05 ratio=0.999962 (target 1.000000 +- 0.02) 
g5_knp_k1280 rc=0 Picard iterations mean=1.827451e+00 | PASS T3 pulse: d(sigma^2)/dt=5.295065e-04 2D=5.208333e-04 ratio=1.016652 (target 1.000000 +- 0.02) 
g5_knp_k12800 rc=0 Picard iterations mean=1.023529e+00 | PASS T3 pulse: d(sigma^2)/dt=5.208134e-05 2D=5.208333e-05 ratio=0.999962 (target 1.000000 +- 0.02) 
g5_kn03_k12800 rc=0 Picard iterations mean=1.023529e+00 | PASS T3 pulse: d(sigma^2)/dt=5.208134e-05 2D=5.208333e-05 ratio=0.999962 (target 1.000000 +- 0.02) 
g5_sq_k128 rc=0 Picard iterations mean=2.000000e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.780463e-03 2D=5.208333e-03 ratio=0.533849 (target 1.000000 +- 0.02) 
g5_aq_k128 rc=0 Picard iterations mean=2.000000e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.753000e-03 2D=5.208333e-03 ratio=0.528576 (target 1.000000 +- 0.02) 
g5_knp_k128 rc=0 Picard iterations mean=2.000000e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.780117e-03 2D=5.208333e-03 ratio=0.533783 (target 1.000000 +- 0.02) 
g5_aq_k12.8 rc=0 Picard iterations mean=2.227451e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.880908e-02 2D=5.208333e-02 ratio=0.553134 (target 1.000000 +- 0.02) 
g5_kn03_k12.8 rc=0 Picard iterations mean=2.003922e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.941880e-02 2D=5.208333e-02 ratio=0.564841 (target 1.000000 +- 0.02) 
g5_sq_k1280 rc=0 Picard iterations mean=1.827451e+00 | PASS T3 pulse: d(sigma^2)/dt=5.295065e-04 2D=5.208333e-04 ratio=1.016652 (target 1.000000 +- 0.02) 
g5_hr_k128 rc=0 Picard iterations mean=2.000000e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.756269e-03 2D=5.208333e-03 ratio=0.529204 (target 1.000000 +- 0.02) 
g5_aq_k1280 rc=0 Picard iterations mean=1.827451e+00 | PASS T3 pulse: d(sigma^2)/dt=5.294702e-04 2D=5.208333e-04 ratio=1.016583 (target 1.000000 +- 0.02) 
g5_hr_k12.8 rc=0 Picard iterations mean=2.000000e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.949960e-02 2D=5.208333e-02 ratio=0.566392 (target 1.000000 +- 0.02) 
g5_aq_k12800 rc=0 Picard iterations mean=1.023529e+00 | PASS T3 pulse: d(sigma^2)/dt=5.208134e-05 2D=5.208333e-05 ratio=0.999962 (target 1.000000 +- 0.02) 
g5_sq_k12.8 rc=0 Picard iterations mean=2.223529e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.872951e-02 2D=5.208333e-02 ratio=0.551607 (target 1.000000 +- 0.02) 
g5_knp_k12.8 rc=0 Picard iterations mean=2.223529e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.873858e-02 2D=5.208333e-02 ratio=0.551781 (target 1.000000 +- 0.02) 
g1_hr_1e2 rc=0 Picard iterations mean=2.060168e+00 | HOPF: max|dE/E| 1.751e-02 top5 5.118e-03 L1 1.557e-03 (E_top 1.7354 exact 1.7344, tau_top-cell 5.09e-04) 
g1_kn03_1e2 rc=0 Picard iterations mean=2.058801e+00 | HOPF: max|dE/E| 1.745e-02 top5 5.070e-03 L1 9.808e-04 (E_top 1.7353 exact 1.7344, tau_top-cell 5.09e-04) 
g3_hr_c10 rc=0 Picard iterations mean=2.526977e+00 | PASS T6 marshak: L1(E)=0.0140 L1(material)=0.0140 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/arms_1010/repo/tests_m1/runs_3a/t6_ref_sn.txt 
g3_kn03_c10 rc=0 Picard iterations mean=2.526397e+00 | PASS T6 marshak: L1(E)=0.0138 L1(material)=0.0138 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/arms_1010/repo/tests_m1/runs_3a/t6_ref_sn.txt 
g3_hr_c1 rc=0 Picard iterations mean=2.006111e+00 | PASS T6 marshak: L1(E)=0.0138 L1(material)=0.0138 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/arms_1010/repo/tests_m1/runs_3a/t6_ref_sn.txt 
g3_kn03_c1 rc=0 Picard iterations mean=2.006111e+00 | PASS T6 marshak: L1(E)=0.0138 L1(material)=0.0138 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/arms_1010/repo/tests_m1/runs_3a/t6_ref_sn.txt 
```

## RESULTS/beams.txt (raw)

```
shd3b hr rc=0 fatal=0 394.84 s wall | Picard iterations mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=2.062000e+01 max=5.100000e+01 | breakdowns=0.000000e+00 |  | solved E <= e_floor cell-solves=0.000000e+00
shd3b sq rc=0 fatal=0 593.63 s wall | Picard iterations mean=4.730000e+00 max=1.400000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.766173e+01 max=3.700000e+01 | breakdowns=0.000000e+00 | plm positivity fallbacks=3.730000e+02 | solved E <= e_floor cell-solves=0.000000e+00
shd3b sqf rc=0 fatal=0 633.51 s wall | Picard iterations mean=5.420000e+00 max=3.300000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.734871e+01 max=3.900000e+01 | breakdowns=0.000000e+00 | plm positivity fallbacks=4.780000e+02 | solved E <= e_floor cell-solves=5.470000e+02
shd3b sqv rc=0 fatal=0 593.32 s wall | Picard iterations mean=4.730000e+00 max=1.400000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.766173e+01 max=3.700000e+01 | breakdowns=0.000000e+00 | plm positivity fallbacks=3.730000e+02 | solved E <= e_floor cell-solves=0.000000e+00
shd3b sqa rc=0 fatal=0 632.84 s wall | Picard iterations mean=5.060000e+00 max=1.300000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.729051e+01 max=3.300000e+01 | breakdowns=0.000000e+00 | plm positivity fallbacks=4.060000e+02 | solved E <= e_floor cell-solves=0.000000e+00
shd3b aq rc=0 fatal=0 458.37 s wall | Picard iterations mean=4.520000e+00 max=1.400000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.729204e+01 max=5.100000e+01 | breakdowns=0.000000e+00 | plm positivity fallbacks=3.600000e+02 | solved E <= e_floor cell-solves=1.600000e+01
shd3b kn01 rc=0 fatal=0 385.58 s wall | Picard iterations mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.825667e+01 max=2.000000e+02 | breakdowns=0.000000e+00 |  | solved E <= e_floor cell-solves=0.000000e+00
shd3b kn03 rc=0 fatal=0 385.34 s wall | Picard iterations mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.824667e+01 max=2.000000e+02 | breakdowns=0.000000e+00 |  | solved E <= e_floor cell-solves=0.000000e+00
shd3b kn1 rc=0 fatal=0 378.13 s wall | Picard iterations mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.828667e+01 max=2.000000e+02 | breakdowns=0.000000e+00 |  | solved E <= e_floor cell-solves=0.000000e+00
shd3b knp rc=0 fatal=0 390.46 s wall | Picard iterations mean=3.460000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.515607e+01 max=2.000000e+02 | breakdowns=0.000000e+00 | plm positivity fallbacks=2.690000e+02 | solved E <= e_floor cell-solves=9.200000e+01
cyl hr rc=0 fatal=0 158.97 s wall | Picard iterations mean=5.000000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.678400e+01 max=1.620000e+02 | breakdowns=3.000000e+00 |  | solved E <= e_floor cell-solves=0.000000e+00
cyl sq rc=0 fatal=0 169.12 s wall | Picard iterations mean=9.790000e+00 max=7.400000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=6.781410e+00 max=9.500000e+01 | breakdowns=3.000000e+00 | plm positivity fallbacks=4.690000e+02 | solved E <= e_floor cell-solves=0.000000e+00
cyl sqf rc=1 fatal=1 198.06 s wall |  |  |  |  | 
cyl sqv rc=0 fatal=0 168.10 s wall | Picard iterations mean=9.790000e+00 max=7.400000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=6.781410e+00 max=9.500000e+01 | breakdowns=3.000000e+00 | plm positivity fallbacks=4.690000e+02 | solved E <= e_floor cell-solves=0.000000e+00
cyl sqa rc=0 fatal=0 198.49 s wall | Picard iterations mean=1.154000e+01 max=1.730000e+02 NON-CONVERGED=0.000000e+00 | inner iterations mean=6.889081e+00 max=1.900000e+01 | breakdowns=0.000000e+00 | plm positivity fallbacks=3.940000e+02 | solved E <= e_floor cell-solves=0.000000e+00
cyl aq rc=0 fatal=0 351.98 s wall | Picard iterations mean=3.380000e+01 max=4.000000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=6.026331e+00 max=1.620000e+02 | breakdowns=3.000000e+00 | plm positivity fallbacks=0.000000e+00 | solved E <= e_floor cell-solves=0.000000e+00
cyl kn01 rc=0 fatal=0 400.22 s wall | Picard iterations mean=4.660000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=8.885622e+01 max=2.000000e+02 | breakdowns=4.000000e+00 |  | solved E <= e_floor cell-solves=0.000000e+00
cyl kn03 rc=0 fatal=0 420.14 s wall | Picard iterations mean=4.870000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=8.867351e+01 max=2.000000e+02 | breakdowns=5.000000e+00 |  | solved E <= e_floor cell-solves=0.000000e+00
cyl kn1 rc=0 fatal=0 413.18 s wall | Picard iterations mean=4.880000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=8.677664e+01 max=2.000000e+02 | breakdowns=7.000000e+00 |  | solved E <= e_floor cell-solves=0.000000e+00
cyl knp rc=0 fatal=0 421.11 s wall | Picard iterations mean=3.399000e+01 max=4.000000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=8.354516e+00 max=2.000000e+02 | breakdowns=0.000000e+00 | plm positivity fallbacks=0.000000e+00 | solved E <= e_floor cell-solves=0.000000e+00
xb20 hr rc=0 fatal=0 423.02 s wall | Picard iterations mean=4.000000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=3.749750e+01 max=8.900000e+01 | breakdowns=0.000000e+00 |  | solved E <= e_floor cell-solves=0.000000e+00
xb20 sq rc=0 fatal=0 1337.52 s wall | Picard iterations mean=6.517000e+01 max=1.770000e+02 NON-CONVERGED=0.000000e+00 | inner iterations mean=7.291392e+00 max=3.000000e+01 | breakdowns=0.000000e+00 | plm positivity fallbacks=6.512000e+03 | solved E <= e_floor cell-solves=2.638000e+03
xb20 sqf rc=1 fatal=1 37.43 s wall |  |  |  |  | 
xb20 sqv rc=0 fatal=0 1263.88 s wall | Picard iterations mean=6.099000e+01 max=1.520000e+02 NON-CONVERGED=0.000000e+00 | inner iterations mean=7.290867e+00 max=3.100000e+01 | breakdowns=0.000000e+00 | plm positivity fallbacks=6.094000e+03 | solved E <= e_floor cell-solves=2.547000e+03
xb20 sqa rc=1 fatal=1 188.65 s wall |  |  |  |  | 
xb20 aq rc=0 fatal=0 406.11 s wall | Picard iterations mean=7.790000e+00 max=1.200000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.457895e+01 max=6.200000e+01 | breakdowns=0.000000e+00 | plm positivity fallbacks=7.360000e+02 | solved E <= e_floor cell-solves=2.140000e+02
xb20 kn01 rc=0 fatal=0 1751.65 s wall | Picard iterations mean=6.950000e+00 max=1.000000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.400849e+02 max=2.000000e+02 | breakdowns=0.000000e+00 |  | solved E <= e_floor cell-solves=1.060000e+02
xb20 kn03 rc=0 fatal=0 1807.17 s wall | Picard iterations mean=7.050000e+00 max=1.100000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.426440e+02 max=2.000000e+02 | breakdowns=0.000000e+00 |  | solved E <= e_floor cell-solves=4.300000e+01
xb20 kn1 rc=0 fatal=0 1817.26 s wall | Picard iterations mean=7.130000e+00 max=1.100000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.424895e+02 max=2.000000e+02 | breakdowns=0.000000e+00 |  | solved E <= e_floor cell-solves=6.300000e+01
xb20 knp rc=0 fatal=0 912.73 s wall | Picard iterations mean=7.970000e+00 max=1.400000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=5.055082e+01 max=2.000000e+02 | breakdowns=0.000000e+00 | plm positivity fallbacks=7.540000e+02 | solved E <= e_floor cell-solves=1.760000e+02
ba0 hr rc=0 fatal=0 485.84 s wall | Picard iterations mean=4.710000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=3.987261e+01 max=1.050000e+02 | breakdowns=0.000000e+00 |  | solved E <= e_floor cell-solves=6.090000e+02
ba0 sq rc=1 fatal=1 85.58 s wall |  |  |  |  | 
ba0 sqf rc=1 fatal=1 45.46 s wall |  |  |  |  | 
ba0 sqv rc=1 fatal=1 81.65 s wall |  |  |  |  | 
ba0 sqa rc=1 fatal=1 66.52 s wall |  |  |  |  | 
ba0 aq rc=0 fatal=0 705.98 s wall | Picard iterations mean=1.544000e+01 max=2.400000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.677720e+01 max=9.500000e+01 | breakdowns=0.000000e+00 | plm positivity fallbacks=1.525000e+03 | solved E <= e_floor cell-solves=1.180000e+03
ba0 kn01 rc=0 fatal=0 778.85 s wall | Picard iterations mean=5.540000e+00 max=7.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=6.632491e+01 max=2.000000e+02 | breakdowns=0.000000e+00 |  | solved E <= e_floor cell-solves=2.439000e+03
ba0 kn03 rc=0 fatal=0 785.17 s wall | Picard iterations mean=5.680000e+00 max=6.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=6.553345e+01 max=2.000000e+02 | breakdowns=0.000000e+00 |  | solved E <= e_floor cell-solves=3.033000e+03
ba0 kn1 rc=0 fatal=0 824.07 s wall | Picard iterations mean=5.700000e+00 max=8.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=6.915088e+01 max=2.000000e+02 | breakdowns=0.000000e+00 |  | solved E <= e_floor cell-solves=3.524000e+03
ba0 knp rc=0 fatal=0 701.39 s wall | Picard iterations mean=1.434000e+01 max=2.000000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.832985e+01 max=2.000000e+02 | breakdowns=0.000000e+00 | plm positivity fallbacks=1.416000e+03 | solved E <= e_floor cell-solves=1.443000e+03
ba20 hr rc=0 fatal=0 555.92 s wall | Picard iterations mean=5.010000e+00 max=6.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=4.549501e+01 max=2.000000e+02 | breakdowns=0.000000e+00 |  | solved E <= e_floor cell-solves=3.900000e+02
ba20 sq rc=1 fatal=1 1719.40 s wall |  |  |  |  | 
ba20 sqf rc=1 fatal=1 95.36 s wall |  |  |  |  | 
ba20 sqv rc=1 fatal=1 1495.94 s wall |  |  |  |  | 
ba20 sqa rc=1 fatal=1 243.33 s wall |  |  |  |  | 
ba20 aq rc=0 fatal=0 952.55 s wall | Picard iterations mean=2.241000e+01 max=2.900000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.666711e+01 max=8.800000e+01 | breakdowns=0.000000e+00 | plm positivity fallbacks=2.230000e+03 | solved E <= e_floor cell-solves=6.910000e+02
ba20 kn01 rc=0 fatal=0 614.07 s wall | Picard iterations mean=4.750000e+00 max=1.300000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=5.580211e+01 max=2.000000e+02 | breakdowns=0.000000e+00 |  | solved E <= e_floor cell-solves=4.160000e+02
ba20 kn03 rc=0 fatal=0 616.93 s wall | Picard iterations mean=4.950000e+00 max=1.300000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=5.337374e+01 max=2.000000e+02 | breakdowns=0.000000e+00 |  | solved E <= e_floor cell-solves=3.760000e+02
ba20 kn1 rc=0 fatal=0 613.23 s wall | Picard iterations mean=4.810000e+00 max=1.300000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=5.500208e+01 max=2.000000e+02 | breakdowns=0.000000e+00 |  | solved E <= e_floor cell-solves=3.330000e+02
ba20 knp rc=0 fatal=0 987.77 s wall | Picard iterations mean=2.062000e+01 max=2.900000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.946460e+01 max=2.000000e+02 | breakdowns=0.000000e+00 | plm positivity fallbacks=2.047000e+03 | solved E <= e_floor cell-solves=6.740000e+02
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/cyl/hr t=2.421e-01
  x in 0.3-0.95: <|E/J-1|> 0.377  max|E/J-1| 0.571  <|F|/|cH|> 0.673  median dangle 4.6 deg  p90 11.0  <f_M1> 0.502 <f_ex> 0.462
  cell flux (derived, clipped): <|F|/|cH|> 0.673  median dangle 4.6  face max |F|/cE 0.65  frac(|F_face| > cE) 0.000
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/cyl/sq t=2.421e-01
  x in 0.3-0.95: <|E/J-1|> 0.379  max|E/J-1| 0.866  <|F|/|cH|> 0.681  median dangle 4.6 deg  p90 12.3  <f_M1> 0.503 <f_ex> 0.462
  cell flux (derived, clipped): <|F|/|cH|> 0.681  median dangle 4.6  face max |F|/cE 0.65  frac(|F_face| > cE) 0.000
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/cyl/aq t=2.421e-01
  x in 0.3-0.95: <|E/J-1|> 0.379  max|E/J-1| 0.646  <|F|/|cH|> 0.673  median dangle 4.6 deg  p90 11.0  <f_M1> 0.502 <f_ex> 0.462
  cell flux (derived, clipped): <|F|/|cH|> 0.673  median dangle 4.6  face max |F|/cE 0.64  frac(|F_face| > cE) 0.000
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/cyl/kn03 t=2.421e-01
  x in 0.3-0.95: <|E/J-1|> 0.377  max|E/J-1| 0.570  <|F|/|cH|> 0.673  median dangle 4.6 deg  p90 11.0  <f_M1> 0.502 <f_ex> 0.462
  cell flux (derived, clipped): <|F|/|cH|> 0.673  median dangle 4.6  face max |F|/cE 0.65  frac(|F_face| > cE) 0.000
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/cyl/knp t=2.421e-01
  x in 0.3-0.95: <|E/J-1|> 0.378  max|E/J-1| 0.629  <|F|/|cH|> 0.673  median dangle 4.6 deg  p90 11.0  <f_M1> 0.502 <f_ex> 0.462
  cell flux (derived, clipped): <|F|/|cH|> 0.673  median dangle 4.6  face max |F|/cE 0.64  frac(|F_face| > cE) 0.000
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/xb20/hr t=2.421e-01  [xb20, x in 0.45..0.95, beam = J > 0.05 Jmax: 2106 cells]
  beam: E/J median 2.930  sum E/sum J 2.821  |F|/|cH| median 2.508  flux-weighted <|F|>/<|cH|> 2.426
  beam: direction error J-weighted mean 4.8 deg  median 3.8  p90 12.6   face |F|/cE > 1: frac 0.0000 (beam) 0.0000 (all)  max 0.94
  dark (J < 1e-3 Jmax, 14108 cells): <E>/Jmax 0.1008  max E/Jmax 1.0909
  x=0.70 shape: L1 |E/sumE - J/sumJ| 0.811
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 4.158  E-centroid 1.0000 (exact 1.0000)  E outside the J>5% band / total 0.396
    peaks y: E [0.965 1.035]  J [1.004]   E max/J max 2.047  width(>50% max) E 0.188 J 0.109
  x=0.90 shape: L1 |E/sumE - J/sumJ| 0.763
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 4.177  E-centroid 1.0000 (exact 1.0000)  E outside the J>5% band / total 0.365
    peaks y: E [0.699 0.879 1.004 1.121 1.301]  J [0.887 0.910 0.926 1.004 1.098 1.113]   E max/J max 2.698  width(>50% max) E 0.391 J 0.281
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/xb20/sq t=2.421e-01  [xb20, x in 0.45..0.95, beam = J > 0.05 Jmax: 2106 cells]
  beam: E/J median 3.072  sum E/sum J 2.904  |F|/|cH| median 2.600  flux-weighted <|F|>/<|cH|> 2.491
  beam: direction error J-weighted mean 4.7 deg  median 3.7  p90 12.7   face |F|/cE > 1: frac 0.0000 (beam) 0.0000 (all)  max 0.94
  dark (J < 1e-3 Jmax, 14108 cells): <E>/Jmax 0.0803  max E/Jmax 1.0982
  x=0.70 shape: L1 |E/sumE - J/sumJ| 0.714
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 3.807  E-centroid 1.0003 (exact 1.0000)  E outside the J>5% band / total 0.331
    peaks y: E [0.965 1.035]  J [1.004]   E max/J max 2.100  width(>50% max) E 0.188 J 0.109
  x=0.90 shape: L1 |E/sumE - J/sumJ| 0.632
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 3.754  E-centroid 1.0010 (exact 1.0000)  E outside the J>5% band / total 0.284
    peaks y: E [0.684 0.887 0.996 1.105]  J [0.887 0.910 0.926 1.004 1.098 1.113]   E max/J max 3.067  width(>50% max) E 0.375 J 0.281
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/xb20/aq t=2.421e-01  [xb20, x in 0.45..0.95, beam = J > 0.05 Jmax: 2106 cells]
  beam: E/J median 2.993  sum E/sum J 2.833  |F|/|cH| median 2.537  flux-weighted <|F|>/<|cH|> 2.428
  beam: direction error J-weighted mean 4.7 deg  median 3.8  p90 12.5   face |F|/cE > 1: frac 0.0000 (beam) 0.0000 (all)  max 0.95
  dark (J < 1e-3 Jmax, 14108 cells): <E>/Jmax 0.1008  max E/Jmax 1.0650
  x=0.70 shape: L1 |E/sumE - J/sumJ| 0.813
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 4.169  E-centroid 1.0006 (exact 1.0000)  E outside the J>5% band / total 0.396
    peaks y: E [0.965 1.035]  J [1.004]   E max/J max 2.073  width(>50% max) E 0.188 J 0.109
  x=0.90 shape: L1 |E/sumE - J/sumJ| 0.765
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 4.141  E-centroid 1.0009 (exact 1.0000)  E outside the J>5% band / total 0.362
    peaks y: E [0.895 1.004 1.098]  J [0.887 0.910 0.926 1.004 1.098 1.113]   E max/J max 2.998  width(>50% max) E 0.375 J 0.281
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/xb20/kn03 t=2.421e-01  [xb20, x in 0.45..0.95, beam = J > 0.05 Jmax: 2106 cells]
  beam: E/J median 2.924  sum E/sum J 2.819  |F|/|cH| median 2.507  flux-weighted <|F|>/<|cH|> 2.424
  beam: direction error J-weighted mean 4.8 deg  median 3.8  p90 12.6   face |F|/cE > 1: frac 0.0000 (beam) 0.0000 (all)  max 0.94
  dark (J < 1e-3 Jmax, 14108 cells): <E>/Jmax 0.1011  max E/Jmax 1.0879
  x=0.70 shape: L1 |E/sumE - J/sumJ| 0.813
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 4.162  E-centroid 1.0000 (exact 1.0000)  E outside the J>5% band / total 0.397
    peaks y: E [0.965 1.035]  J [1.004]   E max/J max 2.044  width(>50% max) E 0.188 J 0.109
  x=0.90 shape: L1 |E/sumE - J/sumJ| 0.765
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 4.180  E-centroid 1.0000 (exact 1.0000)  E outside the J>5% band / total 0.366
    peaks y: E [0.699 0.879 1.004 1.121 1.301]  J [0.887 0.910 0.926 1.004 1.098 1.113]   E max/J max 2.696  width(>50% max) E 0.391 J 0.281
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/xb20/knp t=2.421e-01  [xb20, x in 0.45..0.95, beam = J > 0.05 Jmax: 2106 cells]
  beam: E/J median 3.022  sum E/sum J 2.825  |F|/|cH| median 2.617  flux-weighted <|F|>/<|cH|> 2.427
  beam: direction error J-weighted mean 4.7 deg  median 3.8  p90 12.6   face |F|/cE > 1: frac 0.0000 (beam) 0.0000 (all)  max 0.94
  dark (J < 1e-3 Jmax, 14108 cells): <E>/Jmax 0.1009  max E/Jmax 1.0663
  x=0.70 shape: L1 |E/sumE - J/sumJ| 0.812
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 4.160  E-centroid 1.0001 (exact 1.0000)  E outside the J>5% band / total 0.397
    peaks y: E [0.957 1.043]  J [1.004]   E max/J max 2.149  width(>50% max) E 0.172 J 0.109
  x=0.90 shape: L1 |E/sumE - J/sumJ| 0.762
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 4.163  E-centroid 0.9995 (exact 1.0000)  E outside the J>5% band / total 0.365
    peaks y: E [0.887 0.996 1.113]  J [0.887 0.910 0.926 1.004 1.098 1.113]   E max/J max 4.116  width(>50% max) E 0.180 J 0.281
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/ba0/hr t=2.421e-01  [ba0, x in 0.45..0.95, beam = J > 0.05 Jmax: 730 cells]
  beam: E/J median 7.073  sum E/sum J 5.599  |F|/|cH| median 7.574  flux-weighted <|F|>/<|cH|> 5.862
  beam: direction error J-weighted mean 5.1 deg  median 1.0  p90 4.9   face |F|/cE > 1: frac 0.0096 (beam) 0.0018 (all)  max 78.88
  dark (J < 1e-3 Jmax, 15608 cells): <E>/Jmax 0.1141  max E/Jmax 6.6123
  x=0.70 shape: L1 |E/sumE - J/sumJ| 0.942
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 11.140  E-centroid 0.9806 (exact 0.9732)  E outside the J>5% band / total 0.452
    peaks y: E [0.965]  J [0.980 0.996]   E max/J max 5.894  width(>50% max) E 0.070 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.040
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 12.732  E-centroid 0.9883 (exact 0.9633)  E outside the J>5% band / total 0.516
    peaks y: E [0.957]  J [0.988]   E max/J max 5.715  width(>50% max) E 0.078 J 0.062
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/ba0/sq t=0.000e+00  [ba0, x in 0.45..0.95, beam = J > 0.05 Jmax: 730 cells]
  beam: E/J median 0.000  sum E/sum J 0.000  |F|/|cH| median 0.000  flux-weighted <|F|>/<|cH|> 0.000
  beam: direction error J-weighted mean 2.8 deg  median 3.6  p90 7.6   face |F|/cE > 1: frac 0.0000 (beam) 0.0000 (all)  max 0.00
  dark (J < 1e-3 Jmax, 15608 cells): <E>/Jmax 0.0000  max E/Jmax 0.0000
  x=0.70 shape: L1 |E/sumE - J/sumJ| 1.922
  x=0.70 profile: overlap sum min(E,J)/sum J 0.000  L1 |E-J|/sum J 1.000  E-centroid 1.0000 (exact 0.9732)  E outside the J>5% band / total 0.965
    peaks y: E []  J [0.980 0.996]   E max/J max 0.000  width(>50% max) E 2.000 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.876
  x=0.90 profile: overlap sum min(E,J)/sum J 0.000  L1 |E-J|/sum J 1.000  E-centroid 1.0000 (exact 0.9633)  E outside the J>5% band / total 0.945
    peaks y: E []  J [0.988]   E max/J max 0.000  width(>50% max) E 2.000 J 0.062
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/ba0/aq t=2.421e-01  [ba0, x in 0.45..0.95, beam = J > 0.05 Jmax: 730 cells]
  beam: E/J median 7.054  sum E/sum J 5.598  |F|/|cH| median 7.535  flux-weighted <|F|>/<|cH|> 5.852
  beam: direction error J-weighted mean 5.1 deg  median 1.0  p90 4.9   face |F|/cE > 1: frac 0.0096 (beam) 0.0018 (all)  max 78.88
  dark (J < 1e-3 Jmax, 15608 cells): <E>/Jmax 0.1148  max E/Jmax 6.6820
  x=0.70 shape: L1 |E/sumE - J/sumJ| 0.944
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 11.135  E-centroid 0.9804 (exact 0.9732)  E outside the J>5% band / total 0.452
    peaks y: E [0.965]  J [0.980 0.996]   E max/J max 5.891  width(>50% max) E 0.070 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.040
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 12.722  E-centroid 0.9871 (exact 0.9633)  E outside the J>5% band / total 0.516
    peaks y: E [0.957]  J [0.988]   E max/J max 5.703  width(>50% max) E 0.078 J 0.062
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/ba0/kn03 t=2.421e-01  [ba0, x in 0.45..0.95, beam = J > 0.05 Jmax: 730 cells]
  beam: E/J median 7.073  sum E/sum J 5.603  |F|/|cH| median 7.575  flux-weighted <|F|>/<|cH|> 5.865
  beam: direction error J-weighted mean 5.0 deg  median 1.0  p90 4.9   face |F|/cE > 1: frac 0.0096 (beam) 0.0018 (all)  max 78.88
  dark (J < 1e-3 Jmax, 15608 cells): <E>/Jmax 0.1164  max E/Jmax 6.9674
  x=0.70 shape: L1 |E/sumE - J/sumJ| 0.949
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 11.249  E-centroid 0.9818 (exact 0.9732)  E outside the J>5% band / total 0.456
    peaks y: E [0.965]  J [0.980 0.996]   E max/J max 5.894  width(>50% max) E 0.070 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.044
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 12.829  E-centroid 0.9899 (exact 0.9633)  E outside the J>5% band / total 0.519
    peaks y: E [0.957]  J [0.988]   E max/J max 5.718  width(>50% max) E 0.078 J 0.062
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/ba0/knp t=2.421e-01  [ba0, x in 0.45..0.95, beam = J > 0.05 Jmax: 730 cells]
  beam: E/J median 7.051  sum E/sum J 5.603  |F|/|cH| median 7.533  flux-weighted <|F|>/<|cH|> 5.857
  beam: direction error J-weighted mean 5.1 deg  median 1.0  p90 4.9   face |F|/cE > 1: frac 0.0096 (beam) 0.0018 (all)  max 78.88
  dark (J < 1e-3 Jmax, 15608 cells): <E>/Jmax 0.1177  max E/Jmax 7.1047
  x=0.70 shape: L1 |E/sumE - J/sumJ| 0.954
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 11.298  E-centroid 0.9820 (exact 0.9732)  E outside the J>5% band / total 0.459
    peaks y: E [0.965]  J [0.980 0.996]   E max/J max 5.886  width(>50% max) E 0.070 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.047
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 12.875  E-centroid 0.9890 (exact 0.9633)  E outside the J>5% band / total 0.520
    peaks y: E [0.957]  J [0.988]   E max/J max 5.707  width(>50% max) E 0.078 J 0.062
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/ba20/hr t=2.421e-01  [ba20, x in 0.45..0.95, beam = J > 0.05 Jmax: 797 cells]
  beam: E/J median 2.425  sum E/sum J 2.335  |F|/|cH| median 2.145  flux-weighted <|F|>/<|cH|> 1.999
  beam: direction error J-weighted mean 12.2 deg  median 7.0  p90 31.3   face |F|/cE > 1: frac 0.0113 (beam) 0.0015 (all)  max 78.23
  dark (J < 1e-3 Jmax, 15535 cells): <E>/Jmax 0.0784  max E/Jmax 3.0565
  x=0.70 shape: L1 |E/sumE - J/sumJ| 1.024
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 6.069  E-centroid 0.8618 (exact 0.8959)  E outside the J>5% band / total 0.489
    peaks y: E [0.895]  J [0.918]   E max/J max 2.855  width(>50% max) E 0.078 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.256
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 7.469  E-centroid 0.9191 (exact 0.9590)  E outside the J>5% band / total 0.592
    peaks y: E [0.918]  J [0.980]   E max/J max 2.468  width(>50% max) E 0.156 J 0.070
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/ba20/sq t=0.000e+00  [ba20, x in 0.45..0.95, beam = J > 0.05 Jmax: 797 cells]
  beam: E/J median 0.000  sum E/sum J 0.000  |F|/|cH| median 0.000  flux-weighted <|F|>/<|cH|> 0.000
  beam: direction error J-weighted mean 18.4 deg  median 17.3  p90 22.1   face |F|/cE > 1: frac 0.0000 (beam) 0.0000 (all)  max 0.00
  dark (J < 1e-3 Jmax, 15535 cells): <E>/Jmax 0.0000  max E/Jmax 0.0000
  x=0.70 shape: L1 |E/sumE - J/sumJ| 1.914
  x=0.70 profile: overlap sum min(E,J)/sum J 0.000  L1 |E-J|/sum J 1.000  E-centroid 1.0000 (exact 0.8959)  E outside the J>5% band / total 0.957
    peaks y: E []  J [0.918]   E max/J max 0.000  width(>50% max) E 2.000 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.868
  x=0.90 profile: overlap sum min(E,J)/sum J 0.000  L1 |E-J|/sum J 1.000  E-centroid 1.0000 (exact 0.9590)  E outside the J>5% band / total 0.938
    peaks y: E []  J [0.980]   E max/J max 0.000  width(>50% max) E 2.000 J 0.070
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/ba20/aq t=2.421e-01  [ba20, x in 0.45..0.95, beam = J > 0.05 Jmax: 797 cells]
  beam: E/J median 2.422  sum E/sum J 2.336  |F|/|cH| median 2.152  flux-weighted <|F|>/<|cH|> 1.998
  beam: direction error J-weighted mean 12.3 deg  median 7.1  p90 31.8   face |F|/cE > 1: frac 0.0113 (beam) 0.0015 (all)  max 78.23
  dark (J < 1e-3 Jmax, 15535 cells): <E>/Jmax 0.0779  max E/Jmax 3.0630
  x=0.70 shape: L1 |E/sumE - J/sumJ| 1.018
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 6.026  E-centroid 0.8604 (exact 0.8959)  E outside the J>5% band / total 0.486
    peaks y: E [0.895]  J [0.918]   E max/J max 2.866  width(>50% max) E 0.078 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.248
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 7.376  E-centroid 0.9188 (exact 0.9590)  E outside the J>5% band / total 0.587
    peaks y: E [0.918]  J [0.980]   E max/J max 2.454  width(>50% max) E 0.156 J 0.070
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/ba20/kn03 t=2.421e-01  [ba20, x in 0.45..0.95, beam = J > 0.05 Jmax: 797 cells]
  beam: E/J median 2.422  sum E/sum J 2.335  |F|/|cH| median 2.142  flux-weighted <|F|>/<|cH|> 2.000
  beam: direction error J-weighted mean 12.4 deg  median 6.9  p90 31.7   face |F|/cE > 1: frac 0.0113 (beam) 0.0015 (all)  max 78.24
  dark (J < 1e-3 Jmax, 15535 cells): <E>/Jmax 0.0778  max E/Jmax 3.0066
  x=0.70 shape: L1 |E/sumE - J/sumJ| 1.024
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 6.065  E-centroid 0.8615 (exact 0.8959)  E outside the J>5% band / total 0.489
    peaks y: E [0.895]  J [0.918]   E max/J max 2.852  width(>50% max) E 0.078 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.257
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 7.470  E-centroid 0.9190 (exact 0.9590)  E outside the J>5% band / total 0.593
    peaks y: E [0.918]  J [0.980]   E max/J max 2.467  width(>50% max) E 0.156 J 0.070
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/ba20/knp t=2.421e-01  [ba20, x in 0.45..0.95, beam = J > 0.05 Jmax: 797 cells]
  beam: E/J median 2.436  sum E/sum J 2.351  |F|/|cH| median 2.160  flux-weighted <|F|>/<|cH|> 2.008
  beam: direction error J-weighted mean 12.4 deg  median 7.4  p90 31.6   face |F|/cE > 1: frac 0.0113 (beam) 0.0015 (all)  max 78.24
  dark (J < 1e-3 Jmax, 15535 cells): <E>/Jmax 0.0760  max E/Jmax 2.9969
  x=0.70 shape: L1 |E/sumE - J/sumJ| 1.009
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 5.979  E-centroid 0.8587 (exact 0.8959)  E outside the J>5% band / total 0.481
    peaks y: E [0.895]  J [0.918]   E max/J max 2.872  width(>50% max) E 0.078 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.240
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 7.367  E-centroid 0.9157 (exact 0.9590)  E outside the J>5% band / total 0.579
    peaks y: E [0.918]  J [0.980]   E max/J max 2.550  width(>50% max) E 0.156 J 0.070
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/shd3b/hr shd3b.m1_vsc.00001.bin minJfs 2.600e-13
  x=0.65  E depth 0.000 edge 0.232  Jfs depth 0.000 edge 0.205  ex depth 0.000 edge 0.184  Jfs/ex lit 0.626 E/ex lit 0.593
  x=0.85  E depth 0.000 edge 0.346  Jfs depth 0.000 edge 0.328  ex depth 0.000 edge 0.295  Jfs/ex lit 0.528 E/ex lit 0.608
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/shd3b/sq shd3b.m1_vsc.00001.bin minJfs 2.600e-13
  x=0.65  E depth 0.000 edge 0.204  Jfs depth 0.000 edge 0.205  ex depth 0.000 edge 0.184  Jfs/ex lit 0.626 E/ex lit 0.556
  x=0.85  E depth 0.000 edge 0.344  Jfs depth 0.000 edge 0.328  ex depth 0.000 edge 0.295  Jfs/ex lit 0.528 E/ex lit 0.537
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/shd3b/aq shd3b.m1_vsc.00001.bin minJfs 2.600e-13
  x=0.65  E depth 0.000 edge 0.205  Jfs depth 0.000 edge 0.205  ex depth 0.000 edge 0.184  Jfs/ex lit 0.626 E/ex lit 0.585
  x=0.85  E depth 0.000 edge 0.321  Jfs depth 0.000 edge 0.328  ex depth 0.000 edge 0.295  Jfs/ex lit 0.528 E/ex lit 0.284
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/shd3b/kn03 shd3b.m1_vsc.00001.bin minJfs 2.600e-13
  x=0.65  E depth 0.000 edge 0.232  Jfs depth 0.000 edge 0.205  ex depth 0.000 edge 0.184  Jfs/ex lit 0.626 E/ex lit 0.593
  x=0.85  E depth 0.000 edge 0.347  Jfs depth 0.000 edge 0.328  ex depth 0.000 edge 0.295  Jfs/ex lit 0.528 E/ex lit 0.609
/resnick/groups/carnegie_poc/jingze/arms_1010/run/beams/shd3b/knp shd3b.m1_vsc.00001.bin minJfs 2.600e-13
  x=0.65  E depth 0.000 edge 0.208  Jfs depth 0.000 edge 0.205  ex depth 0.000 edge 0.184  Jfs/ex lit 0.626 E/ex lit 0.631
  x=0.85  E depth 0.000 edge 0.303  Jfs depth 0.000 edge 0.328  ex depth 0.000 edge 0.295  Jfs/ex lit 0.528 E/ex lit 0.518
```
