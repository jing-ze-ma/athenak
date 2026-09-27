---
name: red-giant-top-cell-negative-density
description: "09-11 12:40 ROOT OF THE RT -nan AND (likely) THE 1e16 K VACUUM CELLS: the TOP ACTIVE CELL (i=ie=481) of the open-top red giant gets a NEGATIVE CONSERVED DENSITY mid-stage (u0 before floors); under rt_use_cons PresTempFromEint gives p=NaN, kap=0 at the gate -> dtau_top = 0*NaN -> I_dn NaN down the whole column; the RESCUE (de/e=-0.999) was an ACCIDENTAL ENERGY SINK that let I8 survive. RT guard implemented (uncommitted); with it the same cell breaks the IMPLICIT CONDUCTION solve at 6.1265e5 (131 fallback cells, de=-nan, dt 0.25). Real cure = FOFC on the cubed sphere (step 3) or a top-cell/vacuum treatment"
metadata:
  type: project
---
Evidence: instrumented build (regress/b_rtnan_dbg, -DRT_NANDBG) printed NANDBG TOP m=0 k=2 j=2:
rho_ie=-6.97e-20, T_ie=3.16e10 (ceiling), kap_ie=0, pb=nan, dtau=nan; python confirms
(1-exp(-dtau)) == 0 exactly at dtau>37 so 0*NaN survives to tau 75-90. Sites: two_stream_rt.hpp
:1461 grey, :1681 ck, :1333/:2459 generic (dtau = kap*pb*1e6/EffGravAt).
Fix in the tree (UNCOMMITTED, +144/-25): RTBadState/RTTopDtau helpers, rhoN clamps non-positive
u0 density to dfloor, cells whose (p,T) the EOS cannot form get T_g=0 -> kc=Bb=0 (inert) in all
three kernels, counter + one-time warning. Gate: 150-cycle GV3b bit-identical. Reproduction from
the V9f rst with the I8 input: control reproduces the -nan at 6.8106e5; fixed run has 0 RT NaN
but at cycle 20045 t=6.1265e5 rad_implicit_x1 (conduction.cpp fallback :1096-1110) reports 131
fallback cells with de=-nan at i=481, T* 2e10-3.16e10, rho~1e-17, dt -> 0.2498.
Runs: regress/run_rgctl, run_rgnan, run_rgfix, run_rgfix2; blocks in _analysis_0911/rt_nan_blocks.txt.
Connections: the 1e16 K dt-dip cells of I8 (rho_old<=0, v=-0, energy growing at v=0) are the same
class ([[red-giant-dfloor-energy-diagnosis-wrong]]); FOFC exists for exactly this (trial update
negative density -> first-order fluxes) but is fatal-guarded on the cubed sphere until step 3 of
[[fofc-compatibility-plan]]. Options for the user: (a) same clamp in the implicit conduction (rho
<= 0 -> inert) as a stopgap; (b) FOFC step 3 then enable fofc in the red giant; (c) understand why
the top cell's mass goes negative (outflow flux > content within a stage: WB perturbation
reconstruction at the top? ghost seeding? ceiling momentum?).
Agent cost warning: this chase took 483k tokens / 2319 tool calls (instrumented build + 4 runs).
