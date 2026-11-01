# TASK 2026-11-01 viper -> Caltech (USER): our short bitwise gates move to your GPUs too (viper and Raven queues are blocked)

Use the GPUs freed per NOTE 95ef87ab6. These are short runs, ~150-200 cycles each. Viper and Raven keep their copies queued; the first result wins.

## Builds (H200, fmad off, your usual recipe)
| tag | commit | what |
|---|---|---|
| BASE | fin51-r-fix 65658289c | reference |
| POL | fork/rkgpolar-test-1101 @ 33a4d519f | the USER-approved fix: the RKG super-stage now averages the pole-edge E_r over phi, as the explicit path already does. On sp it changes RKG resistive answers (div B 1.65e-3 -> 3e-11 on sp_test). cs and Cartesian must be untouched. |
| CSV | fork/csvertex-defaults-1101 @ 6029a2aa7 | `<mesh>/cs_seam_fc_nolimit` default ON on fresh starts (recipe::NewDefaults); restarts without the key keep the old behaviour |
| TRC | fork/tracers-defaults-1101 @ 8ccae6636 | the tracers + defaults + ck restart fix + floor mass inherits x_e (you passed T0-T2 on ac98652f) |
| HALO | fork/rkghalo-test-1101 @ 0594b145d | already in TASK 6e5e226a6 part B |

## Gates
Compare rst bytes after `<par_end>` and dhj.mhd.hst exactly; dhj.user.hst floor columns to <= 1e-12. Use any developed W121 cs MHD restart: the r1m p_f4 rst or the T1 rst at 0.25 rot.

| # | test | pass |
|---|---|---|
| G1 | POL vs BASE, key-free, ~200 cycles: (a) RKG non-lean; (b) the r1m production keys (rkg_lean + rkg_dt_safety 0.5) | identical (cs has no polar boundary) |
| G2 | CSV vs BASE: restart a file whose input lacks cs_seam_fc_nolimit | identical, and the log says the key stayed off ("restart input has no key") |
| G2b | CSV fresh start from t=0, ~100 cycles: no key vs `mesh/cs_seam_fc_nolimit=true` | identical |
| G3 | TRC: (a) old r1m rst without tracer keys vs BASE; (b) fresh start without keys = explicit xe_tracer/floor_tracer/valfven_noinfall_fac=4 | identical in both |
| G4 | HALO | as in 6e5e226a6 B |

Report per gate PASS/FAIL in a NOTE. If all pass, we merge POL first (already approved). CSV and TRC defaults go in after your csvertex A/B and the tag-only floor arm.
