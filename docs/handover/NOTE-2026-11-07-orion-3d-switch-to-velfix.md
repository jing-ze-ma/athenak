# NOTE 2026-11-07 Orion -> viper, Caltech (USER): Orion's 3-D copy stopped; 8 nodes now run the velocity-fix 3-D

**USER decision. `prod_seeded` at his settings (job 227350) was stopped cleanly at 1.098 d, with `output/bsg.final.rst` kept, so it can resume. The 17-node takeover (227360) was cancelled.** Raven's copy is now the only run at his exact settings; viper, please say if that changes anything for you.

**Started on the same 8 nodes: job 227441, full 3-D grid, 896 ranks, seeded, with `vel_mode 1`, error_limit 1e-7 and cfl_rad 1 (dt ~219 s).** Run directory `/orion/ptmp/jinma/athpp_repro3d/vf3d/full_e1e-7`. Everything else is his deck and pgen (GM 19.78), with the same dump cadence. Chain of 24 h links, with a 2 TB guard on `output/`.

**Not established:** 1e-7 is chosen for cost; the column showed 1e-6 crashes at this dt, and 1e-7 with the fix has not been checked against 1e-9 past the wedges' early times. SMR with the fix is new here, and the energy budget isn't closed. Expected ~20 h per star-day at first, scaled from the wedge. The wedges at 1e-9, 1e-8 and 1e-7 keep running.
