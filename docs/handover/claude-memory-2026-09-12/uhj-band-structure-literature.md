---
name: uhj-band-structure-literature
description: "What SPARC/MITgcm actually uses for ultra-hot Jupiters: 11 Kataria bins x 8 k-coefficients, 0.26-324 um, in BOTH Parmentier+2018 (WASP-121b) and Tan+2024. Their hosts are all 5500-6500 K. Nobody has published a UHJ SPARC run with an A-star host, which is where the 0.26 um blue edge would break."
metadata:
  node_type: memory
  type: reference
  modified: 2026-08-22T01:30:00.000Z
---

Checked 2026-08-22 to settle whether the 11-band Kataria grid is right for an ultra-hot
Jupiter. It is what the field uses, and it is what Parmentier himself uses.

## Parmentier et al. 2018, A&A 617, A110 (WASP-121b), arXiv:1805.00096

SPARC/MITgcm, Toon+1989 two-stream (delta-discrete-ordinates for the stellar beam,
two-stream source function for thermal), correlated-k with **eight k-coefficients** per
bin. Verbatim: *"When coupled to the GCM, we use 11 frequency bins that have been carefully
chosen to maximize the accuracy and the speed of the calculation (Kataria et al. 2013)."*
Post-processing at higher resolution uses *"196 frequency bins ranging from 0.26 to
300 um"* -- note the blue edge stays at **0.26 um even in post-processing**.

Also states their opacity database *"do not include atomic species such as Fe and Mg
although it has recently been shown that they can contribute significantly to the opacities
in these ultra hot atmospheres (Lothringer et al. 2018)"*. **Our Exo-FMS premixed table DOES
include Fe, Fe+, SiO, TiO, VO** -- we are better off than SPARC there.

## Tan et al. 2024, MNRAS 528, 1016 -- the most recent non-grey UHJ GCM

SPARC/MITgcm, Marley+1999 two-stream with correlated-k, **11 frequency bins** (Kataria
et al. 2015). Post-processing with PICASO at 622 bins, 0.26-267 um. Grids over
**T_eq = 1800-2600 K** and stellar **T_eff of 5500, 6000, 6500 K only**.

## Tan & Komacek 2019, ApJ 886, 26 -- SEMI-GREY, not correlated-k

Used a semigrey scheme to isolate H2 dissociation/recombination. Not a band-structure
reference.

## What this means for our setup

1. **11 bands x 8 g-points, 0.26-324.68 um is the standard for UHJs.** Our choice matches
   Parmentier and Tan exactly. Keep it.
2. **Every published case has a host of 5500-6500 K**, where only 1.3-2.6 % of the stellar
   flux lies bluer than 0.26 um. Nobody has published a SPARC/MITgcm UHJ run with an A-star
   host, where it is 18 % (10000 K) to 32 % (12300 K).
3. **Our omega = 2.06e-5 (3.53 d) does NOT force an A star.** A T_eq = 2500 K planet around
   a 5500/6000/6500 K host would orbit in 0.38/0.64/0.95 days -- six times faster. So either
   the host is hot, OR the planet is not synchronous and omega is an independent rotation
   rate. **The latter is standard practice**: Tan et al. 2024 vary T_eq and rotation period
   as INDEPENDENT grid parameters.

**Recommendation:** set `problem/ck_star_teff` to something in 5500-6500 K, matching the
literature. Then the 11-band grid is well justified, the blue tail is 1-3 %, and results are
directly comparable to Parmentier+2018 and Tan+2024. Only if a genuinely A-star host is
wanted does the 32-band grid (0.20 um edge, premixed table available, ~2.9x the RT cost)
become necessary -- and that would be beyond published practice.
