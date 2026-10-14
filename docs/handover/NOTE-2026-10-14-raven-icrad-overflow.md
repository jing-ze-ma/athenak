# NOTE Raven -> Caltech (10-14): FIN4 int-overflow bug in the fresh-start radiation IC -- please check your He giant runs NOW

src/rad_m1/rad_m1_icrad.cpp:208 tests `vgd_itw.d.extent_int(0) > 1`. vgd_itw.d is the flat per-RANK fused-twin intensity array;
when it holds > 2^31 doubles (He giant fused twin with 2 blocks/rank, and very likely the 1697-cell LARGE He giant even at
1 block/rank), extent_int overflows, the fresh-start radiation IC skips the twin J correction and prints
`ic_rad: no fused twin intensities`. Viper He giant at 2 blocks/rank: 5 % different at t = 0 from every 1-block/rank run.
Restarts are not affected (fresh start only).
Please grep and report:
1. the 200-Rsun He giant FIN4 fresh start (2-GPU link 1 4489484 run.log) -- affected or not;
2. the He giant LARGE gate runs A and C (job 4515555) -- if the line is there, this may be (part of) the BG photosphere mismatch vs E0
   (E0 ran on Raven with 1 small block/rank, unaffected). Include it in the C2/C3 comparison.
The fix is `extent(0)`; Raven will deliver a patched FIN5 (FIN5 build TASK also pushed: TASK-2026-10-14-caltech-fin5-build.md;
build only, no production switch).
