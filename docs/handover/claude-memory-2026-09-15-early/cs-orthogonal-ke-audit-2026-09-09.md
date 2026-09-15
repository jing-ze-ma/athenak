---
name: cs-orthogonal-ke-audit-2026-09-09
description: AUDIT DONE (Opus, read-only) of every place that forms KE/|v|^2 with an orthogonal sum on the cubed sphere — 2 BUG-LIVE (red_giant wall ghost fill, used by prod11; open ghost fill, used by T11/prod12), several LATENT, proposed helpers; NOT yet applied
metadata:
  type: project
---

Basis conventions (verified with file:line by the agent): w0(IVY,IVZ) are CONTRAVARIANT
v^xi, v^eta on unit non-orthogonal tangents with e_xi.e_eta = cos_cell (coordinates.cpp
774-833); u0(IM2,IM3) are COVARIANT m_i = rho g_ij v^j (m2 = d(v2 + c v3), coordinates.cpp
1026-1040); bcc0 is an ORTHONORMAL triple (|B|^2 = sum B_i^2 is right); face fields b.x2f/x3f are
face-normal projections (n_xi.n_eta = -c). Riemann solvers get an orthonormal state
(gnomonic_kernels.hpp 47-100) -> all rsolvers OK. Correct KE: from momenta
0.5[m1^2+(m2^2+m3^2-2c m2 m3)/(1-c^2)]/d; from primitives 0.5 d (v1^2+v2^2+v3^2+2c v2 v3).
|c| <= 0.4755 at a cube vertex -> error scale 33% + cross term up to 1.33x the angular KE.

BUG-LIVE (fix before prod12; also affects prod11):
1. red_giant.cpp `fill` (wall ghost, ~2744-2754; prod11 uses outer_bc=wall): ghost m_i = d v_i
   and E = e + 0.5 d sum v^2 from contravariant w0 -> ghost angular velocity distorted up to
   33% + cross-mixed, ghost e off by (KE_orth-KE_true). ConToPrim runs AFTER the BC so the
   ghost u0 is authoritative. Tens of % at the transonic top of the 24 vertex columns.
2. red_giant.cpp `fill_open_out` (~2711-2727; open top = T11/prod12 config): same, worst case.
LATENT: fill_open inner (~2673), IC (~1473, seed Mach<=0.05), wall_noflux enthalpy (~2412),
open inner passes (~2173, 2242-2281, 2320-2330), MLT shell ske (~1801, mlt off),
srcterms/ismcooling.cpp:225-232 drag mixes contravariant v into covariant m (mesh.cpp:444-447
guard comment wrong; unused here), outputs/derived_variables.cpp v_moments/sgs (output only),
pgen/solid_body_rot.cpp PrimToCons on cs (test), pgen/cubed_sphere.cpp (legacy, dead).
OK: C2P (GnomonicEquiangleRaiseVel over ghosts, defer_cons_floors on cs), history, newdt
(cs/sin_cell), conduction, prolong_prims (refused on cs), fofc/viscosity/most srcterms (refused
on cs), deep_hot_jupiter_rt (Coriolis, ghosts, drag/sponges scale momenta only — they
THERMALISE the removed KE into e_int, flagged as possibly unintended), two_stream_rt/atm_column
(no velocity use).

Proposed fix (apply next session, then rebuild+pin; the running agent owns build_rg_prod):
add to src/coordinates/gnomonic_kernels.hpp device helpers GnomonicKineticFromMom(d,m1,m2,m3,c),
GnomonicKineticFromVel(d,v1,v2,v3,c), GnomonicLowerXi/Eta(d,v2,v3,c); in red_giant.cpp fill /
fill_open / fill_open_out / IC: u0(IM2)=LowerXi, u0(IM3)=LowerEta, et = e + KineticFromVel with
cc = cs ? cos_cell(m,k,j) : 0 (bit-identical off cs); replace the pgen's private CsKinetic;
conserved-side passes -> KineticFromMom; ismcooling drag -> use u0(IM2) directly;
solid_body_rot -> GnomonicEquiangleLowerMom like cs_test.cpp:996.
Full report: .../subagents/agent-a19a562e6f4f84ce7.jsonl (last assistant message).
Related: [[red-giant-sponge-cs-kinetic-energy-bug]], [[cs-nonorthogonal-audit]].
