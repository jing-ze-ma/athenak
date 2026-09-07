---
name: h2-chemistry-quench
description: The tabulated EOS assumes instantaneous H2 dissociation/recombination; recombination quenches above ~1e-3.6 bar, so the upper-atmosphere H2 front is not physical
metadata:
  type: project
---

**Computed 2026-08-26 from the runs' own cells.** The general EOS returns the EQUILIBRIUM
H2 fraction at every (rho, T). That is excellent deep down and wrong at the top.

Rates: three-body `H + H + H -> H2 + H`, k3 = 5.5e-29/T cm^6/s (Palla, Salpeter & Stahler
1983); reverse by detailed balance through the H2 partition function (theta_rot 85.4 K,
theta_vib 6215 K, D = 4.478 eV). t_rec is a fully dissociated parcel re-forming H2, t_dis a
fully molecular one breaking up. ck_limb rot 69; one rotation = 3.05e5 s.

| log10 p [bar] | T [K] | t_dyn (cell) | t_rec | t_dis |
|---|---|---|---|---|
| -1 | 2683 | 7.9e2 | 1.8e-4 | 2.2e-2 |
| -3 | 3318 | 9.0e1 | 7.0e0 | 7.7e-2 |
| -4 | 3337 | 7.4e1 | 7.2e2 | 7.1e-1 |
| -6 | 4154 | 8.5e1 | 1.1e7 | 5.4e0 |
| -7 | 2432 | 3.9e2 | 4.2e8 | 1.3e5 |

**The split is log10 p ~ -3.6** (t_rec = one cell crossing); t_rec = one rotation at -5.3.
Below it the H2 heat pipe -- the reason the general EOS was chosen -- is safe by 4-7 decades.

**Consequence for [[upper-atm-mottling]]:** the cold terminator plumes sit at 1e-7..1e-8 bar,
T ~ 1600 K, n ~ 3e10 cm^-3, where **t_rec ~ 3e10 s, EIGHT decades longer than the 100 s cell
crossing**. Real gas arriving from the dayside would stay atomic (mu ~ 1.3, not 2.3), release
no latent heat and show no density jump. The sinking flow is real; the sharp cold dense blobs
are the equilibrium assumption manufacturing a front. **Do not quote plume amplitudes as
physical.**

**Asymmetry that matters: dissociation stays fast wherever the gas is hot** (5 s at 1e-6 bar,
4150 K). Only the recombination branch -- night side and terminators -- quenches. The dayside
upper atmosphere is fine.

**Nothing cheap rescues it:** the H- route needs electrons and is slower at x_e ~ 1e-6; ISM-
like grain catalysis gives ~2e6 s (~6 rotations), matters over a campaign but not a cell
crossing, and presumes condensates at 1600 K / 1e-8 bar; k3's factor-of-several uncertainty
cannot close eight decades.

## How much of the OBSERVABLE is actually exposed (measured 2026-08-26)

mu comes straight out of a dump as rho k T/(p m_H); fully molecular = 2.32, atomic = 1.25.
On the 1e-6 bar isobar, equatorial band, terminators = |lon| 90 +-20:

| | T | median mu | mu>1.5 | exposed? |
|---|---|---|---|---|
| ck_limb dayside | 4345 K | 1.25 | 2.5 % | **no, fully dissociated** |
| ck_limb terminators | 2273 K | 1.27 | **30.5 %** | partly |
| ck_limb nightside | 1952 K | 1.59 | 56.1 % | yes |
| ck_grav_size terminators | 3229 K | 1.25 | **0.0 %** | no (but only 16 rot) |
| ck_grav_size nightside | 1921 K | 1.67 | 63.6 % | yes |

* **Dayside limb is safe** -- mu = 1.25 to the 5th percentile, equilibrium and frozen-in agree.
* **Terminators are the exposed part**: in ck_limb 30.5 % of the limb ring is molecular
  (13.4 % fully, median 1747 K); those columns would have H up to 1.85x larger if the gas
  stayed atomic, **+9 % in H averaged round the ring at fixed T**.
* **The TEMPERATURE effect is bigger and cannot be bounded statically**: forming H2 at
  f = 0.3 releases ~5e11 erg/g against c_p ~ 1.4e8 erg/g/K. Do NOT claim it is small
  because the mu effect is.
* **Unaffected either way:** domain sizing, isobar radii, floors, dt -- all follow from the
  pressure field, which is what the code solves.

**Inconsistency to keep in view:** `ck_kappa` looks opacity up on (T, p) from the PREMIXED
equilibrium tables and never sees the local composition, so chemical equilibrium is already
assumed on the radiative side. A kinetic f_H2 in the EOS alone would make the two halves of
the model disagree rather than agree.

## DECIDED 2026-08-26: keep as a caveat, do not act

The user's call, after seeing the exposure table: **carry this as a known caveat for the
write-up; do NOT build a kinetic f_H2.** Do not re-propose the advected-scalar model. The
campaign proceeds unchanged -- domain sizing, isobar radii, floors and dt are all
pressure-field quantities and none of them are affected.

**Cheapest passive check, free:** watch whether [[ck-grav-size-run]]'s terminator stays fully
dissociated at 1e-6 bar as it relaxes. If it does, the quench question never touches the limb
observable in the configuration the campaign will actually use.
