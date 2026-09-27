---
name: dhj-photosphere-diagnostic
description: problem/photosphere_dump (a4d30d43) computes the tau=2/3 emitting surface per band and g-point from the run's OWN correlated-k opacity; the photosphere spans 7040 km between the near- and far-IR and the day-night contrast changes 1.7x with band
metadata:
  type: project
---

**Committed a4d30d43**, 2026-09-02. `problem/photosphere_dump = <file>` in
`deep_hot_jupiter_rt.cpp`. Answers "what is the OBSERVED surface" properly instead of
quoting an isobar.

## Why an isobar is not the answer

The emitting level is where vertical tau = 2/3, and it is **not one surface**: deep in a
window between water bands, high in a band centre, and moving across the terminator.
Measured on `sp_dhj_hyd` at 160 rotations:

    band             level        day-night contrast
    1.32-2.02 um     3.1e-2 bar       621 K     <- sees DEEPEST (the H/K window)
    broadband        1.1e-2 bar       634 K
    20-325 um        2.3e-3 bar      1062 K     <- sees HIGHEST

**7040 km of atmosphere** between the deepest and highest band at the substellar point,
and the contrast changes by **1.7x** purely from which band the instrument works in.
Within a single band the window-to-line-core spread runs 30x to 150000x in pressure.
Broadband dayside 1.38e-2 bar vs nightside 8.69e-3 bar.

**The 0.01 bar isobar is a good AVERAGE** (the broadband level is 1.09e-2 bar, r = 1.247
vs the isobar's 1.248 Rp) **and wrong in every individual band.** Zero columns put their
photosphere above the domain top, which independently says the grid is tall enough.

## How it works, and why it is trustworthy

It does NOT re-derive opacity. It calls the same `ck_kappa(...) + rt_kc` the longwave
sweep calls -- correlated-k line opacity plus the CIA, Rayleigh and H- continuum -- and
starts from the same above-domain hydrostatic column the sweep's top boundary term adds.
Per band it reports the g-weighted (emission-weighted) level plus the deepest/highest
g-point; band index CK_NB is the broadband level, Planck-weighted at the local T.

**Needs one RT evaluation first** -- the per-cell tables are allocated lazily, so
`nlim=0` finds them null. Restart and take a cycle:
`./athena -r <rst> -i <input with photosphere_dump> -d . -t 00:02:00`. `-r` and `-i`
together work: the input file is loaded AFTER the restart and overrides it, which is the
only way to add a parameter that is not already in the restart's embedded input.

`ck_kappa`, `ck_planck_frac` and `ck_planck_bands` are now **templated on the view type**,
because under HIP a host mirror is a different memory space and the device-view signature
will not bind. The CPU build compiles either way -- **a host/device space error cannot be
caught by the serial build.**

## POST-PROCESSING TRAP: the code's x2v is a CENTROID, not a stretched cell centre

`coordinates.cpp:1258-1259` sets the spherical-polar cell centres to the **Mignone (2014)
volume-weighted centroid** built from the stretched faces, not to the stretched cell
centre and not to the average of stretched faces:

    x1v = 0.25*(q^2+1)/((1/3)*(q^2+q+1))*(r_r+r_l),   q = r_l/r_r
    x2v = -((sin tr - tr cos tr) - (sin tl - tl cos tl)) / (cos tr - cos tl)

Getting this wrong misplaces the POLE by **1.35 deg** -- largest where the polar cells
are widest and sin(theta) varies fastest across them, ~0 at the equator, so it hides in
any equator-focused check. Found only because the diagnostic prints the code's own x2v
and it disagreed with my reconstruction. **Validate a reconstructed grid against a
coordinate the CODE printed**, not against its own internal consistency.
Cf. [[radial-grid-stretch]] (the stretch itself), [[validate-the-instrument]].

## The visualisation

Artifact "Photosphere Atlas": https://claude.ai/code/artifact/96d83dc5-6806-458a-a50a-f3b580bafdc3
Interactive globe of T / wind on three photospheres and three isobars, plus the band
ladder. Pipeline in the session scratchpad (`load.py`, `eos.py`, `extract.py`,
`gen_data2.py`); the EOS inverter round-trips to 1.0e-3 dex, against the 0.35 dex error a
naive one gives ([[eos-inversion-nan-trap]], [[eos-table-dump]]).

## SECOND POST-PROCESSING TRAP: latitude is NOT uniformly spaced

The user spotted that the published map "seems too far from axisymmetric with respect to
the equator", and it was a rendering bug, not physics. The viewer sampled the field with

    y = (lat - LAT[0]) / (LAT[1] - LAT[0])          // assumes a UNIFORM lat grid

but `f_stretch_theta = 3.0` makes a polar cell 6.8 deg wide and an equatorial one 0.85.
That formula puts the **equator at index 12.50 instead of 31.50** and squashes the whole
field toward the north pole -- which reads exactly as a large north-south asymmetry.
Fixed with a 4096-entry lookup table inverting the real grid; the equator now lands on
31.500 and +/-45 deg map symmetrically to 6.827 / 56.173.

**The measurement that settled it before touching any code:** split T into parts symmetric
and antisymmetric about the equator. At 0.01 bar the antisymmetric part is 6.5% of the
total structure and the two hemispheric means differ by 0.43 K in 2465 K; the grid itself
is symmetric to 8.7e-13 deg. So the run is symmetric to within eddy noise and the defect
had to be downstream.

Note the RAIL STATISTICS were never wrong -- day-night contrast, jet speed and the colour
range are computed in Python on the raw grid, not through the viewer's sampler. Only the
visual mapping was distorted. **A wrong display can coexist with right numbers on the same
page**; check which of the two a complaint is actually about.
Cf. [[validate-the-instrument]].

## THE SURFACE IS CORRUGATED, and a fixed-point ray trace DIVERGES at the limb

The user asked whether the rendered surface accounts for the angle-dependence of the
radius. It did not at first: the sampled VALUES sat on r_phot(lat,lon), but the globe was
drawn as a constant-radius sphere. The relief is not small --

    broadband   1.2333 - 1.2961 Rp   5930 km   5.0 %
    1.3-2.0 um  1.2266 - 1.2763 Rp   4701 km   4.0 %
    20-325 um   1.2493 - 1.3252 Rp   7169 km   6.0 %
    isobars too: 0.95 % at 0.1 bar, 3.90 % at 0.01, 7.81 % at 0.001

-- 5 % of the radius is 14 px on a 620 px globe, so the silhouette is visibly not a
circle. Now ray-traced against r = R(lat,lon), shaded by the true surface normal (radial
tilted by the surface gradient), which is FORESHORTENING not illumination and so is
identical on both hemispheres.

**The fixed-point iteration r <- R(lat(r), lon(r)) does not converge.** Residual over 4
passes: 1035 -> 386 -> 398 -> 434 km, i.e. it drops once then GROWS. Near the silhouette
t = sqrt(r^2 - q) is small, so a tiny change in r swings the sampled lat/lon wildly and
the map's derivative exceeds 1. In pixels on the published globe:

    inner 90 % of the disc   max 0.018 px      <- exact, the iteration is fine here
    0.90 - 0.99              p99 2.30 px
    outer 1 %                median 1.83, max 2.45 px

Fixed with a HYBRID: fixed point inland (q < 0.7225 rdisp^2), and 14-step BISECTION on t
at the limb, where |p(t)| - R(p(t)) rises monotonically so bisection cannot diverge.
Verified on three surfaces: max 0.0022 / 0.0189 / 0.0021 px. Radius quantisation (uint16)
costs 90 m, 1.5e-5 of the relief -- irrelevant.

**Do not report a rendering as accurate because it looks right.** The limb error was
invisible by eye and only appeared when the residual was measured in pixels, and the
non-convergence only appeared when it was measured PER PASS rather than after a fixed
number of them. Cf. [[validate-the-instrument]].

## STILL NOT INCLUDED: the slant path

tau is integrated VERTICALLY, so each point reports the emitting level directly above it.
A real observation looks along a slant path and reaches tau = 2/3 higher near the limb.
That is the standard definition of a photosphere and the right surface to MAP, but it is
not a synthetic image of the disc. Doing that needs tau integrated along rays through the
3-D opacity field, which the diagnostic does not carry -- it stores only the per-column
surface. Stated on the page rather than glossed.
