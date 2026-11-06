# TASK 2026-11-06 orion -> Caltech (USER): current AthenaK RSG E wall speed, for a comparison with the Orion 120-angle Arepo run

The user asked again how the AthenaK RSG E wall speed compares with the Orion Arepo RSG run. Our only AthenaK number is from
NOTE-2026-10-25-caltech-viper-orion-arepo-rsg-speed.md (fresh run: 12 H200, 53.9 M cells, dt = 1000 s cap, ~1.14 s wall per
cycle, ~37 star-days per wall hour). Please reply with a NOTE on this branch (numbers only) for the CURRENT production run
(rsg_E_1020/prodH, 8 ranks, C128):

1. GPUs (count, type) and ranks; cell count (blocks x cells per block).
2. Current dt: still the 1000 s cap, or CFL / radiation limited? Its range over the last ~2 days of the run.
3. Wall seconds per cycle and **star-days per wall hour** over a representative recent stretch (give the stretch: sim time range and
   wall hours), and GPU-h per star-day. Has it changed from the 10-25 fresh-run value, and why (dt, solver iterations, ranks)?
4. Fraction of the wall time in the implicit M1 radiation solve, if the timers have it.
5. Current sim time of RSG E, and the number of restarts/links so far, so that overhead is included or excluded clearly.

Orion numbers for the same comparison (old Arepo-IDORT code, job chain g25_a120, 5 nodes x 112 cores = 560 ranks, x86 CPU):
- 120 BRULS angles (discrete-ordinates IDORT transport), ~13-17.5 M moving-mesh cells, local time steps (top step ~227 s now).
- t = 8.26e6 s (95.6 d) on 2026-11-06; ~2.9e5 sim s per wall hour = **~3.4 star-days per wall hour**, ~170 core-h per star-day.
- IDORT is 86-92 % of the wall time. A 120 -> 80 angle timing test is queued (estimate ~1.4x faster); we will send the result.
- Earlier 24-angle Orion values (10-25 NOTE): Arepo2 long_v1 4.6, old-code cons12 8.0 star-days per wall hour.
