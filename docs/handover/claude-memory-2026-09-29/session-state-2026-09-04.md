---
name: session-state-2026-09-04
description: "START HERE: the granulation artifact is PUBLISHED with the time player and the ray-marched convection-zone plumes; nothing is pending except the user's visual verdict"
metadata:
  type: project
---

Entry point as of 2026-09-04 (last touched at the end of that session; the user asked to
remember and continue next time -- there was no unfinished work, only optional polish). Supersedes [[session-state-2026-08-25]]. Nothing is running
on the cluster.

**Artifact (the one and only — always republish to this URL):**
`https://claude.ai/code/artifact/10828b80-c43a-40ff-8b2c-c6ce558751bd`
Granulation on the corrugated tau=2/3 surface: 80-frame time player + ray-marched
sub-surface plumes. Page is 13.95 MB of the 16 MB ceiling.

**Source now lives at `/orion/ptmp/jinma/Athenak/sunmovie/viewer/` (31 MB, durable).**
Copied out of the session scratchpad on 2026-09-04 because /tmp does not survive.
Build: `cd` there and `python3 build2.py` (needs head2.html, body2.html, app3.js,
atlas.webp, vol/*.webp, frames/meta.json, vol/vol_meta.json) -> granulation_movie.html.
Also there: extract_series.py, mkatlas.py, png_write.py, extract_volume.py, plume_check.py,
atlas.png, and app2.js (the surface-only predecessor). The 80 per-frame surface PNGs were
NOT copied -- atlas.png/atlas.webp hold the same bytes, and extract_series.py regenerates
them from the dumps.

**Data, all from `/orion/ptmp/jinma/Athenak/sunmovie/movie`** (80 dumps, 25 s cadence):
- Surface: `extract_series.py` -> `frames/f0000..0079.png` + `meta.json`, packed into
  `atlas.png` by `mkatlas.py` (8x10 sprite sheet, height 12-bit via `G &= 0xF0`), then
  converted to **WebP lossless** in build2.py: 7.40 -> 4.96 MB, verified byte-identical.
- Volume: `extract_volume.py` (stride 2 -> 40 frames, 50 s) -> `vol/vol_{t,v}{0,1}.webp`,
  96x96x37 over z = 0..1.135 Mm, T' and v_z each normalised by their OWN per-depth sigma
  over the whole series, clipped +-3 sigma, 8 bit, WebP q80. 5.47 MB for all four.
  Per-depth stats print on extraction (sigma_v ~2.6 km/s inside, sigma_T 250-520 K).

**How the viewer works now** (things a future session would otherwise re-derive):
- One WebGL2 `sampler2DArray` holds all 80 surface frames; the surface vertex shader has
  NO attributes (`gl_VertexID` + `texelFetch`), normals from neighbour texels, frames LERPed.
- Plumes: one R8 `TEXTURE_3D` per field per volume frame (80 textures), front-to-back ray
  march in a **half-resolution FBO** composited with premultiplied alpha (ONE,
  1-SRC_ALPHA). 96 steps, jittered. Alpha is `1-exp(-opacity*3*w*dt)` — opacity per Mm, so
  the look does not change with step length. The march is stopped by the lid itself:
  sample the surface texture at (x,z) and break on the first crossing, which gives correct
  opaque-lid occlusion with no depth texture.
- ONE vertical exaggeration now scales lid AND volume together (default 4x, slider 1..30);
  the camera target follows the block centre. Anything else would make the lid cut through
  the fog.
- Transport: play/pause, scrubber, 0.5-4x, Once/Loop/**Bounce** (default, because frames
  80 and 1 are 2000 s apart and uncorrelated).
- `progVol` is built in a try/catch: if the volume shader will not compile the movie still
  runs and the rail says so.

**Budget arithmetic that constrains any future addition:** page bytes = file bytes x 4/3
(base64), so the real budget is 12 MB of files. PNG only reaches ~85% of raw on these
noisy fields; WebP lossless wins ~33% on the surface, WebP q80 ~4x on the volume.

**Are the plumes real? MEASURED** (`scratchpad/plume_check.py`, dump 20, full resolution).
Yes, but with limits worth quoting rather than overselling:
- downflow area fraction 0.51 at the surface -> 0.40 by 0.8 Mm; <v_z> is -1.7 -> -2.3 km/s
  down against a steady +1.7 up, and the fastest 10% of the AREA carries 39% -> 58% of the
  downward mass flux. So sinking material genuinely concentrates and speeds up with depth.
- BUT vertical coherence is modest: corr(v_z(z), v_z(surface)) falls to +0.31 by 0.79 Mm,
  and only ~40% of the deep downflow area sits under a surface intergranular lane. The
  plumes merge and drift; they are not straight tendrils spanning the box.
- the v_z autocorrelation FWHM is 167-208 km at every depth = only 4-5 voxels at the
  volume's 41.7 km resolution. That, plus drawing BOTH signs, is why the first render read
  as mush; the fix was one-sided ("Sinking" by default) + gradient-lit samples.

**User's verdict on the plume view (2026-09-04): "still not very clear but fine."** Accepted
as shipped; do NOT keep polishing it unprompted. The one untried lever is resolution, not
rendering: 128^2 voxels at 20 frames (100 s) instead of 96^2 at 40, which fits the same
budget but halves the plume cadence. Offered and not taken.

**NOT verified:** there is no browser on orion, so the GLSL and the rendered look were
never seen — `esprima` parse + numeric round-trips are all that was checked. If the user
reports a visual bug, that is the first thing to trust.

Related: [[solar-convection-restart-gravity-bug]], [[solar-convection-scaling-wall]],
[[solar-convection-binary-provenance]], [[solar-convection-general-eos]].
