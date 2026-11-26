# NOTE 2026-11-27zzc orion -> Caltech: near-tangent chord fix (OPEN_BUGS #10), bsg + t9l62 arms and sweep cost

From AUDITFIX (orion). USER-approved fix. Idle test nodes only; never stop or share a node with production.
These are two arms added to the combined star check (27zl / 27zu), on the merged tip.

## The change
On curved meshes the closing segment of a ray took chi and S only at its two ends: the upwind point and the node.
A near-tangent chord spans several lateral cells, so an opaque cell lying between the ends was missed.
As a result, E just outside a thick emitter came out about 10 % low (OPEN_BUGS #10).

The fix:
- A chord that spans more than one lateral cell is integrated through the cells holding its sub-points: 4 per cell spanned, at most 16, with chi and S constant in each.
- This happens only where ln chi varies by more than ln 4 along it. A cheap box test runs first.
- It applies on spherical polar and on the cubed sphere.
- It changes only the closure, and every term is positive.

Orion CPU results (vet-nearpath dc1343c05, on caltech-next-3c681d9d1):
- tiny_bsg_rel5 PASS, 3.720e-11 (unchanged).
- quick 27/27.
- Full suite on the pre-gate version: 64/64. Payloads move only in the sp/cs star tests.
- Clump test (tau 20), E at 0-0.5 cell from the edge:
  - sp: -9.6 % -> -7.2 % (0.5-1 cell: -7.4 % -> -3.8 %);
  - cs: -5.1 % -> -1.7 %.
- 1 vs 2x2 MeshBlocks in tile mode:
  - cs: bitwise;
  - sp: 1x2 bitwise; 2x2 at 1e-14, the same as the base.
- Cost on the coarse clump deck (32 radial cells over 1..10, so every near-tangent chord spans 2-3 cells; CPU): ray sweep time +36 % (sp) and +70 % (cs).

The star number is what decides whether the fix is acceptable as is.

## Binaries
- **NB** = your build of caltech-next-3c681d9d1 (27zu).
- **NN** = the same tree + `git am docs/handover/nearpath-1127/0001-nearpath-dc1343c05.patch`. That is one commit, touching rays/sc_ray.hpp, sc_sph_sweep_body.inc and sc_cs_shell.inc.
- Build: CUDA 12.9 + g++13, star_vet.
- Gate: tiny_bsg_rel5 must stay 3.720e-11. Post the md5 and prov.

## Arms (same node type, layout and ranks for each pair)
| arm | binary | star | start | N |
|---|---|---|---|---|
| NB-bsg / NN-bsg | NB / NN | bsg | fresh, template of the tip | 300 |
| NB-t9 / NN-t9 | NB / NN | t9l62 | fresh, template of the tip | 300 |

If 27zl's X-bsg / X-t9 already ran on caltech-next-3c681d9d1 with the same deck and layout, reuse them as the NB arms.

## Report (NOTE caltech-orion-nearpath; files to caltech_results/NEARPATH/<arm>/)
1. Wall time and s/cycle (median of the last 100 cycles), plus the end-of-run line `<rad_vet> vet_gd: rank 0 fenced seconds: ... sweeps ...` of each arm, so the sweep time change NN/NB can be computed.
2. CLEAN WARNING / FATAL / DT_RETRY counts. Also the `lateral reads clamped` / `band-clamped` lines; they must be 0 in both arms.
3. Physics difference NN vs NB at the last cycle: max and mean |rel diff| of rho, T, |v| and E (all cells and rho > 10 dfloor), and the L_top / L_in rows of the last 50 cycles.

-- AUDITFIX (orion), for the coordinator and MERGE
