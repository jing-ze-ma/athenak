---
name: dhj-composition-maps
description: Composition (H2 dissociation, iron ionization, free electrons) mapped from the run's OWN eos_composition.hpp via a standalone dumper; the dayside is atomic and iron-ionized, the nightside molecular and neutral
metadata:
  type: project
---

**2026-09-02.** The Photosphere Atlas now carries composition. The chemistry is NOT
re-derived in Python: `compo.cpp` (session scratchpad) includes
`src/eos/eos_composition.hpp` directly and tabulates on the same (log rho, log T) grid the
run dumps, then Python samples it bilinearly.

**GATE:** it reproduces the run's own `eostab.txt` to **4.8e-7 in log10(p/rho)** with a
**cell-for-cell identical** sub-71 K NaN mask. Build with
`g++ -O2 -I<athenak>/src -o compo compo.cpp` -- it links nothing from AthenaK, the
composition model is a self-contained host-only header with all members public.
`EOSCompositionState` gives xh2, xhii, xe, mu; per-donor gas and ionized fractions are
recovered from `CondensationFactors` and `kmet/nel`, with `nel = xe*rho/(mu*m_u)`.

## What the physics says (sp_dhj_hyd, 160 rot)

* **The dayside takes its hydrogen apart.** Both columns 81% molecular at 0.5 bar; then
  the dayside inversion (2905 K -> 4400 K) dissociates it -- 0.7% at 1e-3 bar, 7e-8 at
  the top -- while the nightside stays 82% molecular to the top.
* **H2 dissociation front** (the OUTERMOST x(H2)=1/2 crossing -- the profile is not
  monotone, the base is only 35% molecular and H2 peaks near the T minimum):
  day **7.2e-3 bar / 1.259 R_p**, night **8.0e-5 bar / 1.336 R_p**. 5.8% of columns never
  dissociate, **all of them on the night side**.
* **Iron is the conductivity switch.** Dayside Fe goes 0.01% -> 99.98% ionized between
  0.5 bar and 4e-7 bar; nightside 0.0000 at every level. x_e differs ~1000x between
  hemispheres. Hydrogen contributes nothing (13.6 eV is out of reach below 4000 K).
* **Iron condensation is NOT a cloud deck**: 0.21% of cells, all above 2.4e-4 bar, 6% of
  them on the 200 K tfloor. Reported on the page, deliberately not drawn.
* Caveat stated on the page: the EOS is an EQUILIBRIUM, and above ~1e-3.6 bar H2
  recombination is 8 decades slower than the flow -- see [[h2-chemistry-quench]].

## The atlas carries BOTH runs now

A **Run** selector switches the whole page between `sp_dhj_hyd` (160.2 rot) and
`sp_dhj_ctl` (100.2 rot); levels, band ladder, composition profiles and stats are all
per-run. Three fields exist only on the magnetized run and are hidden on the other:

* **log10 |B|** -- four decades on ONE emitting surface, 0.02 to 96 G across tau = 2/3
  broadband, median 11 G, from a field that started as a 3 G dipole.
* **B_r, signed** -- where the poloidal field threads out against where it returns.
* **log10 beta** -- minimum **0.96** on the broadband surface: there are pockets of the
  VISIBLE atmosphere where the field sets the pressure, not the gas.

**The two runs are at different epochs, so the atlas is not a controlled comparison** and
says so; the controlled one is [[sp-hydro-vs-mhd-comparison]]. Payload 3.4 MB.
`gen_data4.py` takes `RUNKEY=hyd|mhd` and the two outputs are merged into `dhj_data4.json`.

Artifact "Photosphere Atlas" https://claude.ai/code/artifact/96d83dc5-6806-458a-a50a-f3b580bafdc3
(also brightened this session: ground OKLCH L 0.156 -> 0.248, panel 0.219 -> 0.290).
