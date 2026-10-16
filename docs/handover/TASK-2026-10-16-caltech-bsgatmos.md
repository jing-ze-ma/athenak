# TASK raven -> Caltech: BSG atmosphere-band diagnostic on the Linhao-IC BSG production (bsgli_1014), analysis only

From the Raven coordinator (worker: athpp audit / BSG ATMOS). USER request, 10-16. Post-processing only: no GPU job, no production
change, nothing is resubmitted. CPU on a login or a small CPU node is enough (it reads every hydro_w + the first m1 dump).

## Why
On Raven and Viper we measured the atmosphere just above the photosphere of every old 3-D BSG run (NOTE "BSG ATMOS"). The new
IC-tool setup shows a self-excited 1.1-1.5 d oscillation that saturates at 5-50 km/s; the old Ma-IC runs show a larger infall/refill
motion. Your Linhao-IC production (bsgli_1014, 57.4 d) is the one we cannot reach. Please run the SAME diagnostic on it.

## Files (docs/handover/bsgatmos-1016/, check MD5SUMS)
- atm_runs.py: the diagnostic (below); gate_atm.py: shell means (imported); bin_convert.py: the AthenaK reader (copy of vis/python,
  shipped so the versions match); meta.py: cells per pressure scale height at R_ph.
- Needs python3 + numpy only.

## Run
```
cd docs/handover/bsgatmos-1016
python3 atm_runs.py <RUNDIR> caltech_bsgli > caltech_bsgli.txt        # RUNDIR holds bin/bsg3d.hydro_w.*.bin and bin/bsg3d.m1.*.bin
python3 meta.py <RUNDIR> 2.6542488e27 <R_ph in Rsun printed by atm_runs.py> > caltech_bsgli_meta.txt
```
(If the dumps are not in RUNDIR/bin, give the directory that holds them as RUNDIR; the script uses RUNDIR/bin when it exists, else RUNDIR.)
It writes caltech_bsgli.json next to the script.

## Definitions (what atm_runs.py computes)
- R_ph from the FIRST dump with m1 data: the innermost radius where the shell-mean E <= 1.1 x 4 F_L/c, F_L = L/(4 pi r^2),
  L = median of 4 pi r^2 <F_r> over the inner half of the grid (Eddington tau = 2/3: E = a T_eff^4 = 4F/c; the 1.1 lets an
  isothermal T_ph atmosphere, which only touches 4F/c, still be found). Override with env RPH_RSUN=<value> if needed.
- Shell means over the wedge: <rho> (solid-angle weighted), <v_r> mass-weighted (sum rho v_r dOmega / sum rho dOmega), km/s.
- Bands in r/R_ph: 0.99-1.10 and 1.01-1.10 (atmosphere), 0.50-0.90 (interior). Per dump: max over the band radii of |<v_r>|,
  rms over the band radii of <v_r>, max over the band of |<rho>/<rho>(dump 0) - 1|; also the signed band-mean <v_r> (infall
  indicator: fraction of dumps with band mean < -1 km/s).
- Printed: the max (and when), the median over t >= 1 d, max/median rms, max/median rho swing, signed min/max, infall fraction.
- FFT period of the shell-mean v_r at 1.00 and 1.02 R_ph over t >= 0.5 d, ONLY if the dump cadence is <= 0.4 d (else it says the
  cadence cannot resolve the ~1.3 d acoustic cut-off period; our old runs all had 0.5 d dumps).

## Report back (NOTE-2026-10-16-caltech-bsgatmos.md)
- caltech_bsgli.txt, caltech_bsgli_meta.txt, caltech_bsgli.json (or their contents);
- the binary / branch / deck of the run, its dump cadence, the IC (he_ic_linhao, profile file), sponge keys (he_sponge_*) and floors;
- the wall time of the analysis. Nothing else.
