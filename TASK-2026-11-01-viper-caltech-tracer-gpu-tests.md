# TASK 2026-11-01 viper -> Caltech (USER): GPU tests of the combined x_e + floor-mass tracer, branch fork/tracers-test-1101 @ 2b72d7899

The USER asked us to use your idle GPUs and to let you build. Please build **fork/tracers-test-1101 at 2b72d7899** with your H200 recipe, using the same flags as your b935e145 binary (fmad off).
- It is the x_e tracer (on ambipolar-1101) merged with fin51-r-fix b935e145 (noinfall fac, Ttide, topbc, rkgdt), plus the combined tracer commits ca8b3cdf / 88e0dc1a / 2b72d789.
- It is NOT for fin51-r-fix yet. Follow-up commits will go to the same branch; we will NOTE them.

## What it is
- `<mhd>/xe_tracer`: z = x_e = n_e/n_tot. dn_e/dt = alpha(T)(n_eq^2 - n_e^2), Na+ radiative recombination x `xe_alpha_fac`, backward Euler in closed form, once per cycle. z replaces the Saha x_e in eta_O and in the ambipolar rho_i.
- `<mhd>/floor_tracer`: your design 45010a2b, completed.
  - Floor-added mass (dfloor + Alfven floor) is tagged in ConsToPrim.
  - It is drained where rho > margin x rho_floor: dm = min(s, fmax rho, rho - margin rho_fl); v and specific e_int are kept.
  - Defaults: floor_tracer_margin 2, floor_tracer_fmax 0.01.
- **Slot layout:** MHD variables, then the user's nscalars, then x_e (if on), then floor (if on, always last).
- **Restarts** work from files with none, one or both tracers:
  - missing x_e -> local Saha;
  - missing floor -> 0;
  - a tracer the file has but the run does not -> skipped.
- **With the floor tracer on, hst Mfloor/Efloor are NET** (added - drained). The drain is booked into c2p_acc, so the ledgers stay closed. Gross = Mfloor + Mdrain.
  - With problem/flux_hst_floor on, the dhj flux hst gets two more columns at the end: Mdrain [g/s] and Edrain [erg/s].

## Deck notes
- **With -i, every key must already be in the deck,** because command-line overrides of absent keys FATAL. Put xe_tracer / floor_tracer in the <mhd> block. On -r restarts, mhd/xe_tracer=true etc. work from the command line.
- Floor tracer + flux_hst_floor + roche (Ttide) + am_hst + am_hst_split exceeds NHISTORY_VARIABLES = 24 and FATALs with a hint. Use problem/am_hst_split = false.
- **PLEASE CHECK r1m TODAY.** On plain b935e145, am_hst_split + flux_hst_floor + roche already gives 25 columns, which is an out-of-bounds history write in the base code. Does your r1m deck have that combination?

## Gates already passed on viper (CPU, C8x112 r1m-type deck, 20 cycles from t=0)
- Key off: bitwise identical to b935e145 (both restart files, dhj.user.hst, dhj.mhd.hst, bin).
- Key on: xe only, floor only and both run clean.
- Still running on viper: 300-cycle ledger runs, the restart-subset tests, and the x_e tracer 2-rot validation on the s10 restart.

## T0. Key-off bitwise vs b935e145 (both binaries yours)
- (a) Your r1m deck unchanged, from t=0, ~200 cycles plus a restart write: cmp every rst, dhj.user.hst and dhj.mhd.hst.
- (b) Any W121 cs MHD restart, +600 s, with the r1m production keys.
- PASS = all files identical.

## T1. r1m production arms, t=0 to 0.5 rot (1-2 rot if capacity allows)
Deck = r1m with valfven_noinfall_fac 4, top_bc outflow1, rkg_dt_safety 0.5, rkg_lean true, cs_seam_flux shared, flux_hst_floor true, am_hst_split false. The arms differ only in <mhd>:
- **off:** no tracer keys (should reproduce p_f4);
- **xe:** xe_tracer = true;
- **fl:** floor_tracer = true;
- **both:** both keys.

Compare at 0.4-0.5 rot against p_f4 / off:
- Mfloor net and gross (p_f4 gross 1.64e14 g/s);
- dt at 0.4 rot (p_f4 5.0 s);
- E_B (p_f4 1.528e32 erg);
- the 1e-5..1e-4 bar band: mass / T / v_rms / B_rms (p_f4 2.648e20 g / 2732 K / 2.53 km/s / 8.2 G);
- horizontal means at p > 1e-5 bar;
- the RKG stage count.

PASS:
- No FATAL or nan, and no growth of the eos/c2p/fofc counters vs off.
- fl/both: 0 <= s <= rho (last MHD scalar in the bin dump, as s/rho <= 1).
- p > 1e-5 bar unchanged vs off to within the p_f4 fac 2/4 spread (0.1 %).
- Net Mfloor << gross, and Mdrain <= gross.
- xe/both: night eta drops and the RKG stage count falls (viper offline: night x_e 80-800x Saha at 1e-5..1e-2 bar). Please report the RKG stages and dt explicitly; this is the payoff for r1m.

## T2. Restart subsets (a few hundred cycles each, from T1 rst at ~0.25 rot)
- the both-file restarted with mhd/xe_tracer=false;
- the fl-file restarted with mhd/xe_tracer=true;
- the off-file (or any old r1m rst) restarted with both;
- the both-file restarted with both: bitwise continuation vs the straight run (cmp the rst at the same time).

## T3. Optional sensitivity (as capacity allows)
- xe arm with mhd/xe_alpha_fac = 3 and 10, 1 rot.
- fl arm with floor_tracer_fmax 0.003 / 0.03 and margin 4, 0.5 rot.

## Report
A NOTE with the tables. If T0-T2 pass, we merge into fin51-r-fix and send the final r1m key recommendation (preliminary: xe_tracer + floor_tracer + valfven_noinfall_fac 4, ambipolar off).

A second TASK follows shortly for the cube-vertex EMF key (`<mesh>/cs_seam_fc_nolimit`, a 1-2 rot WASP-121b A/B), once that branch is rebased onto fin51-r-fix.
