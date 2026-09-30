# NOTE 2026-09-30: lhllc / lhlld radial-face fix (lhllc_x1_phi_min), dhj-only default 1

For Caltech / DeltaAI. Branch `lhllc-x1-phi` has been merged into rt-integration. Run data are on viper in
/viper/ptmp2/jinma/w121prod_0929/x1phi_0930 (RESULTS.md, oe_*.txt, cmp_*.txt) and in oddeven_0930.

## Mechanism
In LHLLC (src/hydro/rsolvers/lhllc_hyd.hpp), the contact pressure's velocity-jump term is scaled by
phi = chi(2-chi), with chi = max|v_n|/max c (Minoshima & Miyoshi 2021). The mass flux keeps the full HLLC form.
In the deep dhj envelope the radial Mach number is about 1e-3, so phi is about 3e-3 on x1 faces. This leaves a
2-cell radial velocity checkerboard undamped: odd-even share 0.91-0.96, 97 % of the radial KE, Mdot
alternating in sign between adjacent shells at every depth. Plain hllc removes the mode, but it
lowers the deep horizontal v rms by 7-10 % (low-Mach over-dissipation on x2/x3). LHLLD
(src/mhd/rsolvers/lhlld_mhd.hpp) has the same phi structure.

## Change
- `<hydro>/lhllc_x1_phi_min`, `<mhd>/lhlld_x1_phi_min` (EOS_Data::x1_phi_min, src/eos/eos.hpp):
  on x1 faces only, phi >= phi_min. For lhlld this is applied as the equivalent floor on chi, so the key-off path
  keeps the same compiled arithmetic and stays bitwise. 1 = full HLLC/HLLD velocity-jump term radially; x2 and x3 keep
  the low-Mach fix. x1 is radial on the cubed sphere (gnomonic_kernels.hpp:153) and on sp.
- GLOBAL code default is 0 (off, bitwise). box_convection / He runs, where x1 is the convective
  direction, are unchanged.
- **deep_hot_jupiter_rt sets the key to 1 when the input does not name it**. The pgen runs on restart too,
  and it records the value with GetOrAdd. It prints `dhj: <hydro>/lhllc_x1_phi_min = 1` at startup; check this line in the log.
- Restarts: a restart from ANY older dhj restart file gets 1, because the value comes from the default and not from an input key.
  This was verified: D1 = w1x rst 00600 with no key prints 1, and its restart is bitwise identical to V1. To switch
  it off on such a restart, use `-r file.rst -i key.athinput`, where key.athinput holds `<hydro>` and
  `lhllc_x1_phi_min = 0`. The command line cannot add a key that the restart file does not already hold. Restarts written by the new
  binary carry the key, so later a command-line `hydro/lhllc_x1_phi_min=0` works.

## Evidence (w1x / w10x rot-300 restarts, apudev 2 GPUs, ~5 rot 1x, ~3 rot 10x)
oe = radial odd-even share (1 = pure 2-cell mode), bands 100-1e4 / 30-100 / 10-30 / 1-10 / 0.1-1 /
0.01-0.1 bar. vh = deep horizontal v rms vs the lhllc reference, bands 100-1e4 / 30-100 / 10-30 bar.

| run | oe by band | Mdot alt | radial KE (hst) | vh vs lhllc | jet 0.002-0.02 bar | dt |
|---|---|---|---|---|---|---|
| 1x lhllc (L) | .96 .93 .84 .91 .91 .40 | 1.0 all | 2.0e32 | ref (108/114/189 m/s) | 8069 rms / 18370 max | 13.70 |
| 1x hllc (H) | .10 .04 .02 .08 .21 .07 | <= .25 | 4.9e30 | -9.4/-10.0/-7.6 % | 0.0 % | 13.76 |
| 1x phi_min 1 (V1) | .09 .04 .03 .09 .20 .08 | <= .25 | 8.8e30 | +6.4/+4.0/+1.3 % | -0.6 % | 13.78 |
| 1x phi_min 0.05 | .29 .11 .09 .52 .73 .33 | up to 1 | 1.8e31 | +3.3/+2.2/+0.7 % | -0.3 % | 13.73 |
| 1x phi_min 0.1 | .20 .09 .05 .37 .54 .21 | up to .5 | 1.4e31 | +3.9/+2.6/+0.8 % | -0.4 % | 13.69 |
| 10x lhllc | .85 .71 .42 .32 .73 .37 | 1.0 (.8) | 1.5e32 | ref (76/76/81 m/s) | 6833 / 16271 | 9.56 |
| 10x phi_min 1 | .09 .06 .03 .02 .14 .10 | <= .33 | 1.8e31 | +7.2/+2.5/-1.2 % | +0.3 % | 9.54 |

- p/rho differs from the reference by <= 5e-3 everywhere, in all arms.
- V1's deep vh (100-1e4 bar) rises steadily compared with lhllc: +1.4, +4.0, +5.5, +6.4 % at rot 301, 302.5, 304 and 304.9.
  10x shows the same pattern (+1.4, +4.8, +7.2 %). This is systematic and not yet at equilibrium. Removing the checkerboard lets
  the deep horizontal flow grow, the opposite sign to hllc's over-dissipation. There is no noise member.
- Partial phi floors (0.05, 0.1) leave the mode at 0.1-10 bar, so only 1 is recommended.

## Gates (key 0 / unset)
- dhj key 0 vs pre-merge lhllc: rst data bitwise.
- dhj unset (now 1) vs V1: bitwise.
- He box 30 cycles CPU (box_cpu): hst and bin bitwise.
- MHD deltaai_pkg smoke, 100 cycles, 3 G, lhlld: key 0 bitwise with the old binary, deterministic on repeat. The key-1 arm ran
  cleanly: after 100 cycles the radial KE is 1.39e31, against 1.43e32 with key 0 (2.04e32 at start); dt unchanged.
- tst hydro and mhd CPU tests pass. The rad dhj tests skip (no Exo-FMS tables). Touched files are style-clean.
- Combined gate on the merge: the gated merge a48c5204 (on 2505ca13) has the same src/tst as the pushed merge (on 494540fd, which
  added docs only). It was run from w1x rst 00600 for about 100 cycles (tlim; note that nlim counts ABSOLUTE cycles):
  - key 0 vs the tip: rst data and hydro.hst bitwise;
  - unset: prints 1 and differs, as it should.
  - user.hst columns 12-14 (Efloor/Mfloor/Efloor_rt) differ at 1e-16 even between runs whose restarts are bitwise identical
    (L0 vs L). These columns are a pre-existing nondeterministic reduction, not a change from this branch.
  - He box 30 cycles CPU, tip vs merge: bitwise.
- He box with the key (user question). From hebox_cfl2_0927 H6 at t = 74000 to 79500, cfl 0.6, key 0 vs 1,
  14 restart snapshots from 76000 (x1phi_0930/hekey/ana):
  - the box has only a weak 2-cell v1 component (oe_x1 ~1e-2 in the FeCZ, 0.5 in the top cell), and key 1 cuts
    it to ~1e-4;
  - key 1 changes v1'/vMLT in the FeCZ from 1.108 to 1.059 (-4 %); r(v1',T') is unchanged at 0.235;
  - KE1 is lower by 7 / 11 / 14 % in successive windows, KEh by 1-3 %;
  - dt and the M1 solver counters are unchanged.
  - Verdict: not needed there, and it mildly damps vertical convection. Hence the dhj-only default.
