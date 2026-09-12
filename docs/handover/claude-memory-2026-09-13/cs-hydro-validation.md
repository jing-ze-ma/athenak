---
name: cs-hydro-validation
description: PURE HYDRO cubed sphere, measured end to end -- mass and energy conserve to MACHINE PRECISION including under a SHOCK, Lz cannot and converges at 3rd order, space is 2.3-2.4 (L1) / 1.9 (Linf) and time is 2nd order. Records the FALSE first-order reading a tangential-only refinement gives. START HERE for "does the cubed sphere pass"
metadata:
  type: project
---

Measured 2026-08-31 on `build_cs` (serial, cs_test), iprob=3 rigid rotation for the smooth
cases and iprob=12 blast for the shock, all pure HYDRO, rk2 + plm + hllc.

## CONSERVATION -- mass and energy are MACHINE PRECISION, shock included

Relative drift, initial value from a `time/nlim=0` run of the SAME input, radial
boundaries **reflect** so exactly zero physical flux crosses them:

```
                              mass        energy      Lz
  smooth, unrefined         1.3e-15     -2.7e-15    -1.27e-04
  smooth, SMR (refined)    -5.4e-14      1.4e-14    -3.94e-06  (t=0.25, not comparable)
  SHOCK (blast, p_in 100)  -5.5e-15      1.2e-15     n/a
```

The shock row was **never measured before** -- every gate that printed conserved sums was
a smooth flow. `CSTestBlastCheck` now prints them (437e9b43).

**Lz is NOT machine precision and CANNOT BE**: the gnomonic momentum equations carry
geometric source terms, so Lz does not telescope (already stated in
cubed_sphere_smr.athinput). The right bar is that its drift CONVERGES, and it does, at
**3rd order** -- faster than the scheme itself: -1.268e-04, -1.573e-05, -1.938e-06 at
nx2 = 16/32/64, orders 3.01 and 3.02. The blast carries Lz = 0 by symmetry, so a shock
test with net Lz is still missing.

## CONVERGENCE -- 2nd order in space AND time, everywhere

```
  nx1/nx2      L1(v)        L1(p)       Linf(v)      Linf(p)
   8/16      3.3954e-04   1.0662e-04   1.6990e-03   5.9954e-04
  16/32      6.5117e-05   2.1909e-05   4.9253e-04   1.9746e-04
  32/64      1.2279e-05   4.3565e-06   1.3222e-04   5.3269e-05
  order        2.38/2.41    2.28/2.33    1.79/1.90    1.60/1.89
```

L1 is slightly BETTER than 2nd order; the MAX norm converges at ~1.9 and rising, so the
order holds **everywhere in the active domain**, not just on average. `CSTestConvErrors`
now reports Linf(v)/Linf(p) (MAX-reduced under MPI, not summed).

TIME: at fixed nx2=32, cfl 0.3 / 0.15 / 0.075 / 0.0375 give L1(v) 8.0254960e-05,
8.0252640e-05, 8.0252250e-05, 8.0252160e-05. Richardson on the differences: order 2.57
then 2.12, i.e. **2nd order**, and the temporal error is ~2e-09 against a spatial error of
8e-05 -- **four to five decades below it** at the CFL these runs use.

## THE TRAP: a tangential-only refinement reads FALSE FIRST ORDER

Refining nx2 = 16/32/64 with **nx1 fixed at 8** gives L1(v) orders 2.08 then **1.00**, and
L1(p) 1.62 then 0.51. That is not a defect: the fixed RADIAL resolution is an error FLOOR
the tangential refinement runs into. Refine every direction together, or the test measures
the floor. See [[validate-the-instrument]].

## THESE TWO CRITERIA ARE NECESSARY BUT NOT SUFFICIENT

Also worth gating, and all of them have caught real bugs here:
- **Static/free-stream preservation** -- cubed_sphere_uniform.athinput (rho=1, p=1, v=0
  forever; momenta ~1e-16, KE ~1e-30). A conservative 2nd-order scheme can still destroy a
  hydrostatic atmosphere; conservation and convergence do not see it.
- **Decomposition invariance**, bitwise, over MeshBlock splits and rank counts.
- **HALO norms**, not just domain norms -- a first-order region that shrinks with h is
  invisible to every domain norm. That is exactly what [[cs-wire-fill-wip]] fixed.
- **Long-time secular drift** -- a 3rd-order-small Lz drift still accumulates over the
  many rotations a production run does.
- At a shock, 2nd order is unattainable by construction; the bar there is stability and
  positivity, and note **FOFC is still a startup FATAL** ([[cs-shocks-through-seams]]).
