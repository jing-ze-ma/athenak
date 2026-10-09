# NOTE 2026-10-08 (viper): accretor core fixes are now in rt-integration

For DeltaAI (owner of fork/accretor-1006/1007).

The two core bug fixes from the accretor branch were cherry-picked into rt-integration in the
10-08 merge (branch merge-1008, gated in /viper/ptmp2/jinma/merge_audit_1008/GATE.md):

| accretor commit | rt-integration commit | what |
|---|---|---|
| b226528c | 5120d065 | isothermal pressure p = d c_s^2 in the sp / cubed-sphere geometric sources (IEN read out of range) |
| 1e8d0ff8 | 171daada | periodic NaN scan reads no energy slot for an isothermal EOS |

Same patches (patch-id equal), so `git cherry` shows them as merged; rebasing or merging accretor
onto rt-integration drops/absorbs them without conflict. The rest of the accretor branch
(ry_per_accretor pgen, Hydro::wb_phimax 8cecb89a, data) stays on the DeltaAI branch.

Open, NOT fixed by these: an isothermal HYDRO run on the CUBED SPHERE still reads/writes IEN out
of range in Coordinates::GnomonicEquiangleRaiseVel (called from Hydro::ConToPrim; Kokkos bounds
check aborts at the first ConToPrim, cons index 4 of 4). Spherical-polar isothermal is clean.
GnomonicEquiangleRaiseVelMHD reads u0(IEN) unconditionally too (not run). Do not run isothermal
on the cubed sphere (hydro or MHD) until that is fixed.

Also in this merge (all opt-in / key-off bitwise): he_star_m1 he_gm_column, he_base_heat,
he_wall_inject, Final() View release, he_ic_eint_from_t; mesh f_stretch_r_p2_*; and an
<mhd>/fofc pgen-level default (true for he_star_m1 and box_convection on a FRESH start, as
<hydro>/fofc; restarts without the key keep off).
