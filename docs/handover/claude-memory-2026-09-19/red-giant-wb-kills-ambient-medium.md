---
name: red-giant-wb-kills-ambient-medium
description: "The extended-domain ambient medium died because of the WELL-BALANCED scheme (wb_option=polytropic), not RT/boundary/opacity: WB on -> NaN in 2 cycles with a ~400 g spurious inward force; WB off -> the 400 K background survives and free-falls at exactly g. Fix: a radial WB cutoff (agent adding <hydro>/wb_rmax)."
metadata:
  type: project
---

Decisive A/B on 2026-09-09, runs in /orion/ptmp/jinma/Athenak/red_giant/ (bgonset config:
background rho=2e-22, T=400 K above r=3.89e12, x1max=4.8e12, 6x8x8x320, no RT clipping,
opac_floor 1e-5, direct RT source):

| run | WB | t=527 s | t=875 s | t=3000 s |
|---|---|---|---|---|
| bgnowb | off | e/rho 3.17e10, v_r(i=290) -6.1e3 | -1.0e4 | -3.5e4, no NaN |
| bgwb1 | on, cache every step | NaN by 704 s, e/rho -> 1.42e7 (T~0), v_r -1.3e5 .. -5e5 | | dead |

v_r = -12.7 cm/s^2 * t in bgnowb is free fall at the local g: the physical answer for an
unsupported cold medium (Mdot = the 4e-11 Msun/yr it was specified to deliver).  With WB
on the inward acceleration is ~400 g and the internal energy is drained to below the EOS
table in two cycles, uniformly in angle, from the top down (it looked like an outer-BC
drain, which is why an accretion-boundary change was tried first: bit-identical, no-op).

Everything tried before this was a red herring FOR THIS FAILURE and should not be
re-tried: rt_de_max off, opac_floor, the equilibrium-relaxing semi-implicit update + Newton,
the direct (cancellation-free) RT source, the accretion outer boundary.  Several of those
are still correct changes on their own terms ([[rt-transparent-cell-cancellation]]), but
none of them moved the death time (cycle 35-41, t=2.6-2.8e4) by more than a few percent.

The polytropic WB reference presumes a hydrostatic stratification; a uniform medium is
70 scale heights from one.  R3_nowb (WB off, star only, 1.1 R) runs fine, so a radial
cutoff `<hydro>/wb_rmax` = the star/background join (3.85e12 here) is the targeted fix.
Related: [[red-giant-deep-onset-is-rcb-pileup]], [[red-giant-analysis-opacity-trap]].
