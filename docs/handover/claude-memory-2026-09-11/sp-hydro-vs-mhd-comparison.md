---
name: sp-hydro-vs-mhd-comparison
description: The controlled sp_dhj_hyd vs sp_dhj_ctl pair at 100 rotations -- the field takes 31% of the zonal kinetic energy, 89% of it from below 20 bar, and shows up at the photosphere as only ~2% in radius
metadata:
  type: project
---

**2026-09-02.** `sp_dhj_hyd` and `sp_dhj_ctl` are a **controlled A/B**: identical input
files except `<hydro>` -> `<mhd>`, hllc -> hlld, and `ohmic_resistivity = eos` with
`max_eta = 1e13`. Diff is 22 lines. Bench `/viper/u2/jinma/ATHENAK/bench/compo_ab`.

Both restarted from their own `rst/dhj.00010.rst` and run ~2 min so the photosphere
diagnostic fires; snapshots land at 100.214 and 100.172 rotations -- **0.042 rotations
apart**, and each bin is at the same instant as its own `photo.txt`.

## The result: a magnetic BRAKE, not a sink

From the runs' own history at the snapshot times:

    zonal KE       8.6076e34 -> 5.9633e34   -30.7%
    meridional KE  1.0760e34 -> 6.2074e33   -42.3%
    radial KE      1.0484e33 -> 1.0742e33    +2.5%
    magnetic E                  4.4224e34   = 1.7x the 2.64e34 deficit it opened

**89% of the deficit is in the deepest ten shells**, half of it below 35.8 bar; below
20 bar the magnetized run carries 60-73% less zonal KE. That is where T > 4000 K ionizes
iron, so it is the only place the field couples. **Above ~6 bar the MHD run is FASTER**
(+20% to +375%) -- the energy is moved and stored, not destroyed.

## What it does to the observable atmosphere

* Composition keeps its structure: x(H2) on a fixed isobar differs by at most 3.5%.
* **Every surface sits ~2% further out** (tau=2/3, the 0.01 bar isobar and the H2 front
  alike, +1.7 to +2.6% of R_p, ~200 km) because the MHD run is warmer.
  **The isobar is the control**: its own pressure matches to -0.0%, so the radius shift
  is the atmosphere, not the sampler.
* H2 front 10-17% deeper in pressure; day-night contrast at tau=2/3 **+18%**.
* Fe+/Fe < 2e-4 in BOTH runs at all three surfaces -- the iron switch is higher up, so
  iron is not a discriminant there.

**I first wrote "the composition sits where hydro puts it to well under a percent" and
the table contradicted it** -- the front pressure moved 184% on a geometric mean over all
columns (17% day / 64% night on the median, the honest number). Cf.
[[measure-impact-before-claiming]]. The geometric mean is dragged by the columns that sit
at the model top; use the MEDIAN per hemisphere.

Visualisation: artifact "Magnetic Brake"
https://claude.ai/code/artifact/da2affaa-32b7-4858-96ff-1528e1505746
Pipeline `gen_ab.py` + `shells.py` in the session scratchpad. See
[[mhd-energy-not-from-dumps]] and [[dhj-composition-maps]].
