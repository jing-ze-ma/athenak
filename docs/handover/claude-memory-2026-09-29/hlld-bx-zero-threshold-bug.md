---
name: hlld-bx-zero-threshold-bug
description: 2026-09-14 HLLD/LHLLD BUG (inherited from Athena++): the test `0.5*bxsq < HLLD_SMALL_NUMBER(1e-4)*ptst` compares Bx^2 to the TOTAL pressure, so for plasma beta > 2e4 the rotational discontinuities (** states) are DROPPED regardless of field strength -> with Alfven Mach <~ 1 UNSTABLE (Balsara vortex V~ 1e-3: energies x100-1000, no NaN, global mode, saturates when |Bx| passes the threshold; resolution/CFL/integrator refuted). Fix = 1e-8 for THAT test only (new constant HLLD_BX_ZERO_TOL): bitwise on 3-D linear waves, blast 1e-7, all vortex cases clean; COMMITTED 15e9f790 on lhlld (hlld, lhlld, hlld_uct; 3-D lwave bitwise, blast 1e-8..1e-6, all 4 vortex blow-ups clean); PORT TO polar-average-perf with the RT fixes. Safe sonic Mach scales as sqrt(tol). RELEVANCE: the dhj deep interior has beta >> 2e4 -> every HLLD dhj run evolved the deep field WITHOUT rotational discontinuities (M_A ~ 100s there, so more dissipative, not unstable) - candidate for the sp vs cs deep-field difference; worth one comparison run.
metadata:
  type: project
---
Site: src/mhd/rsolvers/hlld_mhd.hpp ~234, lhlld_mhd.hpp ~310 (check hlld_uct). The other two HLLD_SMALL_NUMBER uses
compare |rho sd sdm - Bx^2| ~ rho c_f^2 and are fine. See [[lhlld-solver]], [[cs-wb-arm-and-ke-gap]].
