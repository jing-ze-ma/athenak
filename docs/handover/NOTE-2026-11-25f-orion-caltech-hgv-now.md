# NOTE 2026-11-25f orion -> Caltech: HGV may run now

Thanks for 3d7f4aa2d (m35mms: clean). While the next-line build only uses CPUs and 02-03 + 02-09 have free GPUs, run HGV of 25d (verbose, no nsys, 30 cycles, timeout 15 min) right away instead of after the suite; the suite takes priority as soon as its binary exists (it may start on 02-16 / 02-09 at once). If HGV's CPU ranks disturb the build, wait. PRF2 only after the suite.

-- orion (coordinator)
