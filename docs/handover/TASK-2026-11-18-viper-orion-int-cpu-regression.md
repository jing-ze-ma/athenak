# TASK 2026-11-18 viper -> Orion (USER plan): CPU regression of clean-cc-int-viper 6f3aa7d07

Branch **fork/clean-cc-int-viper**, head **6f3aa7d07**. It replaces push 18b (42f94c81b) as the clean-cc candidate. It is
b53e35565 (push 18) plus:
- speed8 245ea4418: clean_halo_lists default on (bitwise; falls back on cs).
- L_bot fix 65f075a38: imposed-flux faces (implicit_bc_x1min = flux) get the face stress c E0/3 + 2F/3 + chi d F instead
  of the inner cell's. Without it cell 0 sees ~0.6 of the stress drop of its neighbour at any resolution, the gas falls to
  the wall and L_bot goes negative (R136, RSG). Gate A unchanged (no gate deck uses a flux face).
- physfix2 3ecd48118: write-back margin, clean_stage_floor default false, fref_wacc booking in the cell write-back, your
  closure fix ce55710e3, rad_m1_beam sph_shell KE from the local density, and an energy-conserving gas floor: a cell held at
  e(tfloor) takes the difference from its own radiation E; only what E cannot pay is created and booked (CLEAN floor line).
  The same draw runs before the hydro con2prim (the cs path had no conserved tfloor before). time2_auto_be stays 16.
- realisability opt-ins (all off by default): clean_slave, clean_slave_local (route a': gas gets only the local force in
  slave cells, the rest is booked in am_hist Wp*/Wl*), clean_slave_energy (route a; fails, do not use).
- 201b1f6b7: host team sizes clamp to the execution space (CleanFvRebase aborted on OpenMP with "team size too large";
  the push-18/physfix2 CPU binaries abort on g3_rhse_sp2/sp3 for this reason).
- 6f3aa7d07: rad_m1/clean_nl_tnorm = energy by default (orders identical on the rwt ladders; R136 -9 % s/cycle).

Viper results (GPU unless noted): gate A vs physfix2 c7g 26/30 same (cs_feed1, cs_feed1_x1, gfeed, g8_rhse_cart_pred
differ); restart + forced retry bitwise (7/7); R136 60 cycles from r136.00095 (cycle 494517): L_bot 2.48e39 (old prod
1.98e39; without the L_bot fix -2.5e41), ledger residual -2.2e-4 L_top T, 0 retries; same on Raven CPU.

**Request:** your CPU regression of 6f3aa7d07 against push 18 (b53e35565), as for push 18: pkg_ref gate A, your stage-floor
boxes (sfl/, sflc/) with the t/dt columns and E(gas + rad), and the cs feed box that FATALed at cycle 0 on 42f94c81b.
Please also explain or confirm the g8_rhse_cart_pred difference vs physfix2. Reply with a NOTE. Do not push to clean-cc;
the user decides.

-- viper
