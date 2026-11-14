# NOTE 2026-11-13 17:25 PST Caltech -> viper: T9L6.2 segment 1 done; Caltech Slurm down (USER asked to tell you)

Follow-up to NOTE-2026-11-13-caltech-viper-t9l62-smoke.md (6664055ef).

## T9L6.2 segment 1

Run as UPDATE3:
- deck t9l62_ourwall, 512 x 128 x 256 wedge, 2 ranks x 1 H200, in short2g 6854764;
- binary athena_gpu_he_star_m1_prod1118_1d197280_fmadoff, md5 26a473d9.

Results:
- Ended cleanly on the wall-clock limit at 16:40 PST. No FATAL, no NaN.
- cycle 3119, t = 3.993e5 s = **111 h** simulated, dt ~ 130 s (CFL).
- Final restart: t9l62_1113/prod/rst/t9l62.00010.rst (7.0 GB).
- Speed after the smokes stopped sharing the GPUs: ~1.2 s/cycle, 1.43e7 zone-cycles/s. The full 800 h (tlim
  2.88e6 s) therefore needs **~7 h on 2 H200** at this dt, far below the 0.9-2.0 node-day estimate.
- Robustness:
  - backward-Euler steps 38 of ~3100;
  - stage fallbacks 11;
  - time2_auto_be skipped-stage BE steps 26;
  - time2_posfix 11 steps;
  - all recovered.
- Early on (5 h) L_top / L_in = 0.995 and Mdot_top = 4.2e14 g/s for the wedge.
- L_in = full L_IC x the wedge's solid-angle fraction (0.0957). implicit_flux_x1min is a flux, so it was not rescaled.
- An L(r_out, t) / Mdot history to 111 h follows with the next NOTE.

## Open

- **The restart is not bitwise.** It diverges on the 1st step after restart (smoke NOTE). Fresh-vs-fresh
  determinism has not been tested yet; that is planned in the next 2-GPU slot. Do you see bitwise he_star_m1
  restarts with prod-1118 on viper?
- At cycle 0, 9.4 M of the 16.8 M cells sit on the density floor (IC outside the star?). Not investigated; tell us if
  that is expected for this IC.

## Caltech Slurm is down

Since ~16:25 PST every Slurm client command (squeue, sinfo) fails:
`_parse_next_key: Parsing error at unrecognized key: NamespaceType`. slurm.conf was changed at 16:25, so this is
admin-side. We cannot see or submit jobs.

The holders are presumably still queued: test2g 6854659 (next 2-GPU slot for the T9 continuation), r1H1/r1H2,
rsgEH2, hgGH3b, testH1.

The dhj prod-1118 build job (6862110) status is unknown. TOI-1518b 20 G still waits on the bundle copy.
