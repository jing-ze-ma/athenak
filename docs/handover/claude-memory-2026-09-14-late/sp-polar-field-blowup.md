---
name: sp-polar-field-blowup
description: THE SPHERICAL-POLAR CONTROL HAS ITS OWN NUMERICAL FIELD BLOW-UP AT THE POLE -- LOCALISED 2026-09-04 (later session) to the m = nx3/2 ODD-EVEN CHECKERBOARD of B_r (and B_phi) in the POLAR ROW, onset rot 0.6 at the INNER radius, growing in a STILL axisymmetric fluid (v ~ cm/s): the FIELD OPERATORS alone drive it. The B_theta polar average (use_polar_average_b) does NOTHING (sp_polavg). The polar ghost exchange is EXACT (1e-16). Arms running: both polar flags, HLLE, no resistivity; plus 200-cycle CPU operator ablations from rst 10
metadata:
  type: project
---

**Found 2026-09-04 while comparing cs_rot32 (matched resolution, corrected rotation) against
sp_dhj_ctl through the dt trough.** The two grids agree on total magnetic energy to 1-6 %
through rot 0.6. Then:

    rot    1-ME (r)    2-ME (th)   3-ME (phi)   total       growth per 0.1 rot
    0.70   1.74e32     3.72e32     3.80e32      9.26e32     1.28x
    0.80   9.67e33     4.19e32     2.48e33      1.26e34     13.57x   <-- 56x in B_r
    0.90   1.28e34     4.40e32     3.06e33      1.63e34     1.29x
    1.10   1.45e34     4.94e32     3.51e33      1.85e34     1.04x    (saturated)

A 56x jump in the RADIAL component in a tenth of a rotation, saturating at once, is not a
jet winding a seed field (that builds B_phi, slowly). It is a numerical instability hitting
its ceiling. **RETRACTED the same hour: "77 % of B_r energy within 5 deg of the poles at the inner
boundary."** That was computed from the DUMP's cell-centred bcc, and the instrument check --
integrate 0.5 B^2 dV over the dump and compare with the history -- FAILED: 0.007 of the
history's radial ME at rot 2 (0.35 at rot 0), while theta came out at 0.89, so it is not a
volume-element error. It is [[mhd-energy-not-from-dumps]]: bcc is the mean of the FIELDS,
the history integrates the mean of the ENERGIES from the faces, and at rot 2 they differ by
140x in B_r. A quantity that captures 0.7 % of the energy cannot localise it. WHERE the
blow-up lives is UNKNOWN until it is measured from FACE fields (a restart file), not a dump.

**What the 140x DOES establish, more robustly than any localisation**: 99.3 % of sp's
radial field energy at rot 2 cancels in a cell average, i.e. B_r alternates sign at the
GRID SCALE. A grid-scale sign-alternating radial field that appears in 0.1 rot and
saturates is numerical, wherever it sits. (The record's "8.8x off in B_r" was this same
ratio earlier in the run; it has grown to 140x.)

**Both cubed-sphere runs (nx=32 AND nx=64) sit at 0.82-0.84 of sp at rot 0.69 -- i.e. where
sp was BEFORE its jump -- and never reproduce it.** The cubed sphere has no polar axis.

## What this contaminates

* [[sp-hydro-vs-mhd-comparison]]: "the field takes 31 % of the zonal KE, 89 % of it from
  BELOW 20 bar" -- a strongly magnetised deep state that appeared in 0.1 rot as grid-scale
  B_r. Treat as suspect until the face-field localisation is done.
* Every "cs vs sp" comparison of MAGNETIC quantities after rot ~0.75, including the dt
  yardstick: sp's dt trough at rot 1.1 (1.20) is deeper than cs's (3.35 at rot 0.64) plausibly
  because sp is carrying a spurious strong grid-scale field. dt is the wrong yardstick.
* The premise that sp is "the clean control". It is clean in HYDRO. In MHD it has a
  polar-axis instability of its own, consistent with [[cs-mhd-instability-characterized]]'s
  finding that the low-beta instability is NOT gnomonic-specific (sp_lowbeta reproduces it).

## What it does NOT change

The cubed sphere's own low-beta seam inconsistency on the omega = 0 reproducer
([[cs-mhd-low-beta-divergent]]) stands. Two grids, two different numerical weak spots: the
seam on cs, the polar axis on sp. Neither is a trustworthy MHD control for the other above
beta ~ 1 near its weak spot.

Open: whether cs's KINETIC energy deficit (~40 % of sp from rot 0.4 on, BEFORE the field
jump) is also polar on the sp side, or is the cs fallback's dissipation. Measured next.

## Instrument lesson, again

The polar localisation was written and retracted within the hour because the check
"does my volume integral reproduce the history" was run AFTER the claim instead of before.
On this code, ANY energy localised from a dump must first reproduce the history total;
if it does not, the dump quantity is not the one the history integrates.


## LOCALISED PROPERLY, from FACE fields (restart dhj.00001.rst, rot 10), 2026-09-04

Restart layout identified rather than assumed (ng=2, nvar=5 reproduce data_size exactly);
block geometry paired from the dump at the same rotation. Summing 0.25(bf^2 + bf^2) -- the
history's own formula -- WITHOUT volumes, so the comparison cannot be a volume error:

    theta-ghost B_r energy density        pole-touching blocks   interior blocks
      theta-LO ghost layer                2.189e9  (4 blocks)    1.211e6  (12 blocks)
      theta-HI ghost layer                2.163e9  (4 blocks)    1.162e6  (12 blocks)
      ALL ACTIVE CELLS, all 16 blocks     2.898e9

**Polar ghosts / interior-block ghosts = 1834x** (per block ~5500x). The eight polar-ghost
layers hold MORE radial field energy density than the entire active domain. It is the
POLES, not theta block boundaries in general -- the interior theta-ghosts are the null
control and they are quiet. So this is a defect of the spherical-polar POLAR BOUNDARY
CONDITION for the face field, leaking into the adjacent active cells, which is what the
history's 13.6x jump (active cells only -- history.cpp loops nmb*nx3*nx2*nx1 from is/js/ks)
records.

Also measured: only 8 % of adjacent radial-face pairs have opposite-sign B_r over the
whole domain, so the "grid-scale alternation" is confined to where the energy is, not
global.

**STILL OPEN, and not to be over-read**: my volume-weighted active-cell integral gives
0.13 of the history's radial ME (theta 1.32, phi 0.85) with spherical volumes, and the
Cartesian-volume hypothesis was REFUTED (0.000; sp_dhj_ctl was built at d25d8792,
2026-09-01, AFTER the 985faa22 volume fix). The 7.7x normalisation gap is unexplained --
plausibly the volume of the first active cell at the pole, which dominates -- and it does
not affect the ghost-vs-ghost contrast above, which uses no volumes.

## MECHANISM (2026-09-04): B_r is MULTI-VALUED at the polar axis, and nothing enforces otherwise

At the axis B_r must be single-valued, i.e. independent of phi. From the rot-10 restart's
face fields, for every pole-touching block, azimuthal mean vs azimuthal RMS of the m!=0 part
of B_r on x1 faces, summed over radius:

    first GHOST row past the pole   |mean_phi| ~ 2.0     rms(m!=0) ~ 4.3e2    ratio ~ 200
    first ACTIVE row at the pole    |mean_phi| ~ 2.3     rms(m!=0) ~ 4.3e2    ratio ~ 190

**The radial field at the axis is ~200x more non-axisymmetric than axisymmetric.** The ghost
rows reproduce the ANTIPODAL block's active rows exactly (block 0 ghost == block 8 active), so
the cross-pole exchange is correct; the defect is in the active cells at the axis.

Why it is permitted: the polar face-field BC (`bfield_bcs.cpp`, `case BoundaryFlag::polar`)
sets ONE value -- `x2f` at the pole face, as the mean of its neighbours -- and touches neither
x1f nor x3f in the theta-ghosts, while every other BC type fills all three through all ng
layers. The routine that would keep B single-valued at the axis, `PolarAzimuthalAverageBxBy`
(`bfield_bcs.cpp:307`, declared in bvals.hpp:198), has been COMMENTED OUT at its only call
site (line 216) since the polar boundary was introduced in 9fb17531 (2026-06-18) -- it has
never run. The EMF average `PolarAzimuthalAverageEr` IS live (mhd_corner_e.cpp:396, 478), but
it constrains E, not B: an m!=0 B_r at the axis, once seeded, is not removed.

This is spherical polar's analogue of the cubed sphere's seam problem: each grid has one
singular structure, and each has an unenforced constraint there. Neither is a clean MHD
control for the other.

**The test that would confirm it** (not run): re-enable `PolarAzimuthalAverageBxBy` on
sp_dhj_ctl's input from scratch and check that (a) the rot-0.75 jump in history ME
disappears, and (b) this m!=0/m=0 ratio at the axis stays O(1). If (a) holds, the
"31 % of zonal KE in the field" result in [[sp-hydro-vs-mhd-comparison]] must be redone.

## THE TEST IS RUNNING (queued 2026-09-04): `bench/sp_polavg`

`<mesh>/use_polar_average_b` (committed, default OFF) re-enables `PolarAzimuthalAverageBxBy`.
`sp_polavg` = sp_dhj_ctl's input + that flag, FROM SCRATCH, restarts every 2.0e4 s, six
chained apudev slots (11386143-48). Same HIP binary as the cs arms plus this flag.

Gates passed on the polar blast (run/pole_blast reduced to 32x32x16, ONE block in r -- the
polar boundary refuses more): flag on, 1 rank == 2 ranks bitwise (the MPI_Allreduce path);
flag changes the answer (not a no-op); axis m!=0/m=0 ratio 0.1 with it on AND off -- the
blast is healthy at the pole either way, so it proves harmlessness, not efficacy. Only the
dhj run shows the disease.

**Read it with `bench/pole_axis_diag.py <restart> <any sp_dhj_ctl dump>`** (block geometry is
static, any dump works) and the history: the questions are (a) does the rot-0.75 ME jump
disappear, (b) does the axis ratio stay O(1). Outcomes: both -> sp_dhj_ctl's MHD after 0.75
was the pole and sp-hydro-vs-mhd must be redone; ratio drops but jump survives -> B_r needs
its own m=0 average; neither -> mechanism is elsewhere.

## A RESISTIVE-SPECIFIC GAP in the polar EMF average (found 2026-09-04)

`MHD::EField` calls `CornerE` first -- which ends with `PolarAzimuthalAverageEr()` on
`efld.x1e` -- and only THEN `AddResistiveEMFs(b0, efld)` into the same array. So the polar
azimuthal average of E_r covers the IDEAL EMF only; the resistive part (eta J_r) is added
afterwards, un-averaged, and at the axis it is free to be multi-valued. Spherical polar
only: on the cubed sphere the resistive EMF goes in BEFORE the seam average. LIVE in
production (`ohmic_resistivity = eos`, explicit, `use_rkg_sts = false`), including in the
sp_polavg test that was already running when this was found.

Fix: `<mesh>/use_polar_average_eresist` (default OFF) re-applies the average after the
resistive add; idempotent on the ideal part. Kept OFF by default so sp_polavg stays
reproducible on its own binary. **The next sp arm should run with BOTH polar flags on.**

The audit that found it, for the record: rotation source and outer-x1 ghost KE are in the
pgen (path-independent); resistivity builds NO cell-centred field of its own (consumes bcc0
from RaiseVelMHD, so the stretched-grid interpolation flows through); the RKG STS sub-loop
calls the same ApplyPhysicalBCs -> BFieldBCs, so the B average would run there too; the
fallback diagnostic counts Riemann faces only, by design.

### Trap: the polar blast CANNOT gate resistivity

run/pole_blast with `ohmic_resistivity = constant` collapses dt from 1.6e-4 at cycle 0 to
4e-9 at cycle 1 -- the explicit diffusion limit dx^2/eta with dx = r sin(theta) dphi -> 0 at
the pole -- then crawls 35 cycles and grinds to NaN (eta 1e-3 AND 1e-5). Any "gate" on it
compares NaN with NaN and reads IDENTICAL. Production survives only because EOS resistivity
is negligible at the pole. Gate resistive polar changes on the dhj input itself, reduced
(nx2=16 nx3=32, one block) and run serially for a few cycles.

### Both polar flags ACT (gated on the right binary, 2026-09-04)

One step from sp_dhj_ctl's rot-10 restart on a binary that contains the flags
(functionally verified: sample pole-face x2f 0.3885 -> 0.2496), data-region bytes
differing from base: `use_polar_average_eresist` 50,599 of 98 M (a localised EMF-at-the-pole
change); `use_polar_average_b` 21.5 M of 98 M (the whole state moves once the pole-face
B_theta is rewritten). The earlier "IDENTICAL" gates were on a binary built BEFORE the
flags existed -- see [[validate-the-instrument]], the binary-vintage trap.

## VERDICT (2026-09-04, end of session): the B_theta average does NOT stop the jump

sp_polavg (from scratch, `use_polar_average_b = true`, eresist OFF), history:

    rot 0.700   totME 9.635e32   1-ME(r) 1.789e32     (control 0.70: 9.26e32 / 1.74e32)
    rot 0.709   totME 1.038e33   1-ME(r) 2.324e32
    rot 0.846   totME 1.595e34   1-ME(r) 1.250e34     (control 0.80: 1.26e34 / 9.67e33)

The jump happens in full, marginally LARGER than the control. The average was verified
acting (21.5 M bytes moved per step on the diseased state), so this is a real negative:
constraining the pole-face B_theta to a single Cartesian vector does not remove the
non-axisymmetric B_r growth at the axis. The mechanism is not sourced through B_theta at
the pole face, or not only.

**Next arm**: (1) both flags on (add `use_polar_average_eresist`, the resistive E_r at the
axis is still un-averaged) -- cheap, same input; (2) if that also jumps, an m=0 constraint
on B_r ITSELF in the polar cells (azimuthally average x1f over the ring adjacent to the
axis each step), which nothing currently enforces; (3) localise the FIRST cycles of the
jump from sp_polavg's own restarts (every 2e4 s, so the 0.70-0.85 window is covered --
unlike sp_dhj_ctl's 10-rot cadence) with `bench/pole_axis_diag.py`, to see whether the
m!=0 B_r appears in the ghost row or the active row first.


## 2026-09-04, LATER SESSION: THE MODE IS LOCALISED (bench/pole_axis_spectrum.py, pole_corner_probe.py, pole_corner_flow.py, pole_map_check.py)

**sp_polavg (use_polar_average_b = true) is IDENTICAL to the control**: 1-ME 2.3e32 -> 1.25e34
between rot 0.71 and 0.85, same as sp_dhj_ctl. The B_theta average is NOT the cure. On its
restarts (every 2e4 s), the active-row m!=0/m=0 ratio of B_r at the axis goes
3.6, 21, 65, 357, 432, 394 (rst 09..14) -- exponential then saturating; the history ME jump
IS the saturation.

**WHAT the mode is.** Full-ring FFT of x1f B_r in the first active polar row:
Nyquist share of the m!=0 power 0.0 % at rot 0.2-0.5, 4 % at rot 0.525 (rst 08), 96 % at
rot 0.656 (rst 10), 99.6 % at rot 0.8. Sudden ONSET, not a linear mode growing from
round-off. Power peaks at i = 0 (inner radial boundary), 60-70 % in the first 8 cells.
In the polar active row B_r AND B_phi are 98-99 % pure Nyquist at EVERY radius, while
rho, M_r, M_th, M_phi, E carry 0.00 of it: a PURE-B odd-even mode, a CT null structure
(alternating-sign radial hairs closed by a tiny B_phi checkerboard). The flow there is
still and axisymmetric through the onset (v_phi ~ 8 cm/s, rms 1-3 cm/s; v_r a uniform
60 m/s outflow; rho rms 1e-6 relative) while B_r rms goes 4e-3 -> 2.2 -> 52 G on a 3 G
mean (rst 04, 08, 10). So the FIELD OPERATORS ALONE drive it: CT corner-E (GS05),
HLLD's dissipation of a transverse-B jump (zero at v_n = 0, B_n = 0 -- a tangential
discontinuity is an exact HLLD steady state), the resistive EMF, the polar E_r average,
the inner-x1 user BC (mirrors B_r with area ratios, REVERSES B_theta/B_phi, ghost div B
uncorrected -- src/pgen/deep_hot_jupiter_rt.cpp:2485).

**What is EXONERATED (measured):** the polar ghost exchange -- B_r(ghost,k) =
+B_r(active, k+N/2), B_theta and B_phi = -(mirror), to 3e-16 on rst 04, both poles.

**Arms launched (all from the sp control, ONE change each):**
bench/sp_polboth (both polar flags, 6 slots), bench/sp_hlle (rsolver = hlle),
bench/sp_noeta (resistivity removed). And bench/sp_cpu_ablate/{ctl,hlle,noeta,nopolavg,
eresist}: 200 cycles on the serial CPU binary from sp_polavg rst 10 (t = 2.0e5 s), where
the mode is pure and still linear; measure the Nyquist amplitude growth with
pole_corner_probe.py. NOTE a command-line override on a RESTART only works for parameters
present in the EMBEDDED input; otherwise pass -i with a modified input (main.cpp: -i
overrides the restart's parameters).

**Read [[cs-rotation-source-bug]] first: cs never shows this because it has no pole, not
because it is healthier.** [[validate-the-instrument]]: the bcc dumps CANNOT see this mode
(cell mean of alternating faces); use the restart's face fields.

## 2026-09-04, END OF SESSION: THE MODE IS THE (-1)^(i+k) CT CHECKERBOARD -- INVISIBLE TO bcc

Signed Nyquist amplitude of x1f B_r per radial FACE i on sp_polavg rst 10 (rot 0.656):
+51.1 -50.4 +49.7 -48.7 ... : it alternates in RADIUS as well as in phi. x3f B_phi is the
same (-1)^(i+k) pattern at +-20 G. The CELL-CENTRED B_r (mean of faces i, i+1) is only
0.3-1 G, i.e. 100x smaller. So the mode is the classic CT null checkerboard that
cell-centred quantities cannot see: the RIEMANN SOLVER (which reconstructs bcc) neither
sees nor damps it. That is why the 200-cycle HLLE ablation (bench/sp_cpu_ablate/hlle,
rst/dhj.00011.rst, t 2.0e5 -> 2.003e5) GREW it by 2.4 % in 300 s -- e-fold 1.25e4 s, the
same rate as the GPU run. Only operators built from FACE fields touch it: the resistive
EMF (CurrentDensity, src/diffusion/current_density.hpp, j2 uses dx1*B_r differences in k
and dx3*B_phi differences in i), the CT corner-E derivative terms, the polar E_r average,
and the inner-x1 user BC. It starts at i = 0 -> the inner boundary is where the pattern
is seeded (the BC mirrors B_r with an area ratio and REVERSES B_phi; an odd-even pattern
in i is exactly what a mirrored ghost can feed).

**WHERE TO PICK UP (in order):**
1. bench/sp_cpu_ablate/{ctl,noeta,nopolavg,eresist}/rst/dhj.00011.rst (were still running,
   ~23 s/cycle on one core). Measure with `python3 ../pole_corner_probe.py <rst>
   ../sp_polavg/bin/dhj.mhd_w_bcc.00000.bin` and compare the B_r Nyq amp at i=2 against
   the start value 4.97e1 (rst 10) and HLLE's 5.09e1. If noeta grows the same -> the
   resistive EMF is irrelevant and the driver is CT/BC; if noeta is flat or damped ->
   the resistive curl at the polar row is anti-diffusive on this mode (check the
   (-1)^(i+k) response of j2/j3 in current_density.hpp by hand).
2. GPU arms sp_polboth / sp_hlle / sp_noeta (6 apudev slots each, jobs 11409414-57):
   verdict at rot 0.7-0.85 in dhj.mhd.hst column 11 (1-ME); control jumps 1.7e32 -> 9.7e33.
   Expect sp_hlle NOT to cure it (the ablation says so).
3. The fix should act on the FACE fields at the source, e.g. a divergence-free ghost fill
   at the inner x1 boundary (the commented-out div_rest block at
   deep_hot_jupiter_rt.cpp:2500) or an odd-even filter of x1f/x3f in the polar row; gate
   on the signed-Nyquist amplitude, not on bcc or history ME.
4. cs production reruns cs_prod_hyd_rot / cs_prod_mhd_rot started on apu1 at ~18:30
   (see [[cs-dhj-production-retry]]).

## 2026-09-05: MECHANISM FOUND -- the GS05 corner-EMF upwind terms grow the mode at ~v_r/dr

**All ablations grow the mode by the SAME factor** (bench/pole_nyq_growth.py, rst 10 -> 11,
300 s, every radius, B_r and B_phi alike -- a pure eigenmode): ctl 1.02411, hlle 1.02408,
nopolavg 1.02411, eresist 1.02411, noeta 1.02476 (resistivity DAMPS it by 0.06 %, so eta is
irrelevant). GPU arms agree: sp_polboth, sp_hlle, sp_noeta all jump to 1-ME ~1e33-1e34 by
rot 0.8-0.9. e-fold = 300/ln(1.02411) = 1.26e4 s.

**Prediction from the GS05 stencil (src/mhd/mhd_corner_e.cpp "emf3" kernel):** the corner
e2 derivative terms e2x3 - e2cc are fed by the PHI-FACE Riemann EMF, which contains v_r x
B_phi(face) and therefore SEES the (-1)^(i+k) face checkerboard that bcc (and hence e2cc and
the x1-face flux) cannot. Working the stencil through for that mode with v_r > 0 gives a
corner e2 = -(v_r C/2)(-1)^(i+k), a DOWNWIND-biased EMF, so dB/dt = +(v_r/dr) B: growth at
v_r/dr. Measured at the inner active cell: dr = 7.95e7 cm (stretch map), v_r = 6.1e3 cm/s
(pole_vr_probe.py) -> dr/v_r = 1.30e4 s vs 1.26e4 s observed. The mode is closed ONLY in the
polar row because the polar ghost row is the k+N/2 image of the row itself, which the
Nyquist mode is invariant under; in interior rows a single-row mode leaks in j.

**Decisive ablation launched:** `<mhd>/bs_emf` (new flag, plain four-face average in place
of GS05 on any grid) as bench/sp_cpu_ablate/bs, paired with ctl2 on the SAME rebuilt binary
(build_dhj_cpu; the earlier ablation binary had a different md5 -- binary-vintage trap).
Prediction: bs growth factor ~1.000. UCT (mhd_corner_e_uct.cpp, MZ21) is DEAD code: not in
CMakeLists, its arrays are not declared -- not a ready fix path.

**Decisive ablation READ (2026-09-05, rst 10 -> 11, same rebuilt binary):** ctl2 1.02411 (digit-
identical to ctl -- the binary is faithful); **bs_emf 0.999** (growth GONE; the GS05 corner
terms ARE the driver). But BS is NOT a fix: its dt collapsed 1.97 -> 0.36 s in 300 s, i.e. the
plain average destabilises something else (it does not reduce to the upwind flux in 1D). The
fix candidate is `<mhd>/polar_emf_diss` (Rusanov face-field dissipation on e2 in the polar
rows, local fast speed) = bench/sp_cpu_ablate/diss; GPU arm staged in bench/sp_diss.

## 2026-09-05 FIX VERIFIED ON CPU, committed 1cabe85c (pushed to fork)

`<mhd>/polar_emf_diss = true` (bench/sp_cpu_ablate/diss): checkerboard 50 G -> 3e-4 G in
300 s (factor 1e-5), dt 1.968 unchanged, mass/tot-E identical, floor counters identical,
ME2/ME3 unchanged, 1-ME down 15 % (the polar hairs' energy). GPU arm **bench/sp_diss**
(jobs 11416367-72, six chained apudev slots, from scratch, sp_dhj_ctl input + ONLY this
flag): gate = 1-ME in dhj.mhd.hst at rot 0.7-1.0 stays ~1e32 (control jumps to 9.7e33 by
0.8) AND the signed-Nyquist probe (bench/pole_nyq_growth.py) on its restarts. If it passes,
redo sp-hydro-vs-mhd with the flag ([[sp-hydro-vs-mhd-comparison]] is suspect), and make
the flag default ON for use_polar_boundary MHD runs.

**GPU ARM PASSES (2026-09-05 ~21:30): bench/sp_diss** at rot 0.771: 1-ME 6.77e31 (control
9.67e33 at 0.80, 140x), dt 1.69 normal, KE3 6.13e33 = control's pre-jump value; signed-Nyquist
probe on rst 12: 1e-6 G (diseased run 50 G). The polar fix WORKS on the GPU from scratch.
Still to do: let it reach rot 1.0+ (slots 11416370-72), then (a) make polar_emf_diss default
ON when use_polar_boundary && MHD, (b) rerun sp MHD long enough to redo
[[sp-hydro-vs-mhd-comparison]] -- the control's KE3 9.9e34 at rot 7 (7x its hydro run) is
almost certainly the blow-up's imprint.
**GATE PASSED (2026-09-05 22:xx): sp_diss rot 1.13, 1-ME 7.6e31 flat (control 1.39e34 at 1.0),
dt 1.21-1.37 = the control's own range, KE3 6.0e33.** Slots exhausted; the run stops at 1.13 rot.
The sp polar-axis blow-up is CLOSED by `<mhd>/polar_emf_diss` (1cabe85c). Not yet default-on.

**2026-09-05 late: DEFAULT ON (3147de4f, pushed) and the sp MHD RERUN LAUNCHED:
bench/sp_mhd_diss** (sp_dhj_ctl's input + explicit polar_emf_diss=true, from scratch, HIP
binary of 3147de4f, apu1 2 GPUs, jobs 11420831-32 chained 24 h; the control did 101 rot in
one slot). Purpose: redo [[sp-hydro-vs-mhd-comparison]] against sp_dhj_hyd (its horizontal
KE 1.9-2.3e34 at rot 6-22). Gate: 1-ME stays ~1e32 past rot 0.8; KE3 vs sp_dhj_hyd.
