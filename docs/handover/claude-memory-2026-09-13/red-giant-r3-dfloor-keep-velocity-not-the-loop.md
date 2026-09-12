---
name: red-giant-r3-dfloor-keep-velocity-not-the-loop
description: "R3 RESULT 09-11 07:20: dfloor_keep_velocity is CORRECT (audited, 3 defects fixed, gates pass, inert on the healthy star) but NOT the loop: R3_bface died bit-for-bit like R2 (cycle 19448, t 5.8667e5); collapse cell rho 3e-18 (3x dfloor, NOT floored), v ~7e11 cm/s = 1e5 x escape (9e6); energy created in RKUpdate at KE-dominated cells (ei = E - KE cancellation); RT column already NaN (T(i-1)=nan from negative ei under rt_use_cons). NEXT: Newtonian <hydro>/vceil (none exists; eos_vceil is SR/GR only) + NaN guard in the rt_use_cons T precompute -> R4"
metadata:
  type: project
---
Agent A audit fixes (working tree): general_hyd.cpp writes cons(IEN) when dfloor_fv<1 on the
non-deferred path (was silent internal-energy creation); ideal_hyd.cpp + general_hyd.cpp skip
the momentum/dfl_fv writes on the only_testfloors (FOFC) pass; eos.cpp fatals if the flag is
set in <mhd>. No separate counter (eos_dfloor counts the same events). pin9 = pin8 + these.
Gates: pin7 vs pin9 (flag off, rt_bface=true) bit-identical; flag on vs off on the healthy star
(8400 cycles) bit-identical. Run dirs GA*/N2_bface/N3_dfl under red_giant/.
R3_bface (196173, pin9, rt_bface+keepv): collapse cell (2,9,7,395) gid 26 r 3.641e12 rho 3.04e-18
T 8.3e6 v=(6.3e11,-4e11,7.3e11); first crossing [after_RKUpdate+conduction] gid 26 (9,7,396)
T 3.16e10 rho 4.8e-16 de/e=1; crossings 158 vs R2 150, same first crossing gid 85 cycle 18567.
R3_scan (196156, rt_bface OFF by mistake) = "bface off" control: collapses from 5.791e5.
I7_dfloor/ prepared, NOT submitted. Escape speed at 3.6e12 (1 Msun) ~ 9e6 cm/s; the runaway
cells move at 1e11-1e12, so a ceiling of a few 1e7 is physically free and kills the
KE-cancellation (KE/ei ~ 1e9 makes ei from u0 garbage). See [[red-giant-runaway-source-rt-stale-w0]].
