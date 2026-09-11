---
name: cs-radial-unification
description: 979edada -- the FULL cs-vs-sp radial audit, and cs now treats the radial direction exactly as sp (volume CENTROID, position-aware x1 reconstruction always, weighted radial B, dhj pgen uses stored x1v/xx1f). What was audited clean, what changed, the gates. READ before any new cs radial claim
metadata:
  type: project
---

**2026-09-06, the user's ask:** "the radial treatment should be the same for both in hydro,
mhd, resistivity and deep_hot_jupiter_rt; cs should use the volume centroid, not the
midpoint." Done in 979edada (after [[cs-stretched-resistive-rcm-bug]] c5c85e3b and
[[cs-stretched-source-term-bug]] 13a97399).

**Audited SHARED already (one curvilinear branch, pcoord areas/volumes/dx):** hydro/mhd
update, mhd_ct, resistivity_ct/update, hydro/mhd newdt, resistive dt, flux_correct_cc/fc,
prolongation/restriction (area and volume weighted), history volume, divb diagnostic.

**What DIFFERED and is now identical to sp:**
1. Coordinates cs kernel used the stretch-mapped MIDPOINT r_c for dx2, dx3, dxface,
   areaedge, r_cm; x1v (reconstruction, pgen) is the CENTROID. Now RadialCentroid(r_l,r_r)
   (grid_stretch.hpp) everywhere, as sp's x1v_.
2. x1 reconstruction: cs took GridPiecewiseLinearX1 only when stretched; now ALWAYS (the
   centroid is off-centre even on a uniform grid). `reconstruct` governs ANGULAR sweeps
   only on cs -> the per-direction split the ppmx arm needed exists now for free.
3. Resistive energy-flux radial weights, cell-centred radial B in ConsToPrim (ideal +
   general) and in GnomonicEquiangleRaiseVelMHD: unconditional on cs.
4. dhj pgen: its 13 cs radial positions take x1v_/xx1f_ like the sp branch.
5. cs_test: cell centres are centroids.

**Gates:** sp dhj production input 3 cycles old vs new: hst + bin dumps IDENTICAL. cs dhj
input: mass 3.46013e26 -> 3.46161e26 (IC at centroids). cs_test iprob=9 uniform: within
10 % of before; STRETCHED 8-cell: L1(B) 4.4e-5 -> 1.19e-5, L1(v) 1.14e-4 -> 6.4e-5 = the
uniform numbers. iprob=11 stretched: Ohmic 0.998, curl residual 2.49e-4 = uniform.

**Consequence:** ALL cs numbers change slightly (centroid), so every cs baseline in memory
predating 979edada is stale by ~1-10 % in test errors and 0.04 % in dhj mass. The
rcmfix arm was resubmitted on this binary (job 11442863; 11442237 cancelled at rot 0.4).
Still by DESIGN different: sp src3 flux form vs cs state form; pwb subtraction sp-only.
