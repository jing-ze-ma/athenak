---
name: sun-convection-viz
description: "The run/sun solar-convection box -- 3D volume renderings and the tau=2/3 granulation surface, published as the \"Buoyancy Box\" artifact. Also the MEASURED verdict that the structures above the convection zone are internal gravity waves and NOT shocks, and the fact that the 512-cubed convection data is GONE"
metadata: 
  node_type: memory
  type: project
  originSessionId: bb0f18f2-0ed4-4a24-964f-a406a6150a1e
  modified: 2026-09-03T23:36:07.393Z
---

**2026-09-03.** Artifact **"Buoyancy Box"**, https://claude.ai/code/artifact/2050a98a-4d63-4202-9796-d43fe4eb5fc7
-- interactive: the tau = 2/3 surface (temperature / vertical velocity, 41 frames), the
volume-rendered interior (41 time frames + 36 azimuths, drag to rotate), and the
waves-vs-shocks evidence.

## THE DATA THAT IS GONE

The **512^3** Cartesian low-Mach convection run (`cs_test`-unrelated; `src/pgen/convection.cpp`,
Leidi+2024) COMPLETED on 2026-02-07 -- job 5618651, `Terminating on time limit`, tlim 4000,
10 h 21 m on 32 nodes -- and **its snapshots were deleted**. `sacct -j 5618651 -o WorkDir`
says `run/convection`, which today holds only the two log files. The 256^3 twin (job 5618658)
wrote to `run/convection_256`, which no longer exists at all. Nothing on viper, raven or the
HPSS archive at `/ghi/r/j/jinma` (which holds one .bash_history and nothing else).

**Resolution is read off the log's FIRST dt**, since the log does not print cell counts:
dt(cycle 0) = 1.171875e-3 = CFL 0.3 x (2/512) / 1.0 EXACTLY. At 256 it would be 2.34e-3.
Root grid "4 x 4 x 4 MeshBlocks" is 64 blocks of 128^3.

The setup survives: `/raven/u/jinma/ATHENAK/athenak/run/convection/convection.athinput`
(256^3 there) and `src/pgen/convection.cpp` (a USER problem, so `-D PROBLEM=convection`).
Cost to regenerate 512^3 to tlim 4000: **331 node-hours**, and 400 dumps of `hydro_w` at
512^3 is **2.1 TB** -- which is presumably why it was deleted. Design the outputs for the
figure, not for the original run. **The user says the dir is on ORION, not here.**

## WAVES, NOT SHOCKS -- measured, on run/sun

The `run/sun` box (`solar_convection`, 128^3, 41 dumps, t = 0..40 ks) has convection below
z = 0.65 Mm (entropy flat, Mach 0.065-0.2) and a stably stratified layer above it. The
striations above the interface are **internal gravity waves**, on four independent measures:

    max compression -div v * dx/cs   0.048 below the sponge (0.30 anywhere) -- a shock
                                     needs ORDER 1, so it is 20x too weak
    |div v| / |curl v|               0.13, i.e. the field is 87-89% SOLENOIDAL
    entropy vs pressure              rms ds / rms(dp/p) = 1.9-4.5 and ANTI-correlated;
                                     sound is adiabatic and carries no entropy
    (dp/p) / (gamma dr/r)            0.29-0.57;  an acoustic wave sits at exactly 1.0

The trap: max Mach DOES exceed 1 late in the run (1.51), but **every supersonic cell is at
z = 1.19-1.29 Mm** -- inside the sponge damping layer, against the top boundary, in the
thinnest gas. Numerics, not physics. Clip the sponge (top 20%, z > 1.05 Mm) out of any figure.
A concentric ring at z ~ 0.95 Mm looks like a front but is a VORTEX ring at a plume head:
local div/curl 0.71 against 0.13 background, but corr(dp,drho) = +0.55 and amplitude ratio
0.288, both far from the +1 and 1.0 an acoustic front requires.

## The tau = 2/3 surface

tau integrated from the top down with the pgen's OWN `get_kapr` (H--like T^9 branch,
bound-free, electron scattering above 1e4 K, 1e-2 floor); T = p/(rho*Rgas) with the run's
Rgas = 1.38e8, validated against the pgen's quoted base state (recovers 13,740 K vs its
stated 13,718 K). Found in all 16,384 columns at a mean height of **0.60 Mm** -- which lands
on the 0.65 Mm interface located independently from the entropy profile.

    T at tau=2/3     3494-4349 K, contrast 2.83% rms; the hot end matches Teff ~ 4160 K
    v_z at tau=2/3   only 41.5% of the AREA is rising: discrete upflow cells inside one
                     connected sinking network, whose lanes match the dark lanes in T
    corrugation      +-20 km against a 10.2 km CELL -- barely two cells

**Colour maps come from `run/plotsun.py`, the user's own script for this run**: `seismic` for
vertical velocity at +-1 km/s, `jet` for temperature. One deliberate deviation: plotsun.py
spans 0-6000 K for a vertical slice, but the SURFACE spans only 3494-4349 K = 14% of that
ramp, so temperature limits are fitted to the surface. The +-1 km/s velocity scale is kept
exactly and saturates 13.7% of the surface.

**Do NOT build the corrugation into an exaggerated 3D landscape.** It needs x45 to be visible
and that dresses up a barely-resolved quantity as structure; use it as relief SHADING under
the colour.

## Renderer

No pyvista / VTK / yt / plotly and no ffmpeg-on-PATH on viper (`module load ffmpeg` works).
The volume rendering is ray-marched emission-absorption written directly in numpy
(`volrender.py`: rays clipped to the box, 210 samples, trilinear via
`scipy.ndimage.map_coordinates`, front-to-back compositing). Two things that made it read:
a transfer function with a DEAD ZONE (below 0.30 of full scale nothing is drawn, and the rest
maps onto the saturated wings, skipping the colour map's pale centre, which otherwise fogs
the box to grey), and an opacity BOOST above the interface, because wave fronts are thin
sheets that a ray crosses in two samples and they otherwise vanish beside the fat plumes.
The field rendered is v_z / rms_h(z) -- density falls four decades, so without that
normalisation the convection is all you can see. For a MOVIE the normaliser must be ONE FIXED
profile, or it divides out the very growth the movie is about (wave rms at 0.9 Mm grows 8x
between 4 and 40 ks).

Rendering is single-threaded numpy on a 256-core login node: split frames across ~12 worker
processes with a stride, and STAGGER their launch -- a simultaneous numpy import storm over
NFS kills several with `ImportError: numpy C-extensions failed`.
