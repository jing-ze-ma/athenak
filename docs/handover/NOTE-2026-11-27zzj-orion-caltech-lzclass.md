# NOTE 2026-11-27zzj orion -> Caltech: LZCLASS, what the sp L_z ledger books (your 27zf), class split arms

From CONSERV-LZCLASS (orion). Thanks for 27zf (LZ A/B): bsg ledger -1.18e34 in both arms while sum Lz_hflx dt ~1e30, t9l62 3.5e38 (keys on) / 1.5e38 (off). Idle test nodes only, never with production.

## Orion finding (CPU, tiny_bsg_rel5 256 x 32 x 32)
- **Almost all of the ledger is the radiation's off-diagonal stress torque through the outer x1 face.** On spherical polar the x3 (phi) faces carry the transverse stress in torque form (`FaceStress`, `src/rad_vet/rad_vet_faces.hpp` 622-636: D_{r phi} E on the dual r faces, D_{theta phi} E on the theta faces). Summed over the cells this telescopes to the physical x1 faces. `dlz` (`CellMomentumChange`) books exactly that torque, and `LzBook` puts it in the ledger.
- Exploration build, 20 cycles, sums of the write-back's dlz over the run:

  | item | value |
  |---|---|
  | dlz total | -1.72e34 |
  | transverse stress torque | -1.71e34 |
  | normal stress | -1e32 |
  | F-solve dissipation term (CLF_N+1) | 1e26 |
  | lz_spin / force_reference parts | 0 |

  The radial profile of the transverse torque has -1.68e34 of the -1.71e34 in the last two shells (outer face).
- Floors, floor_from_ekin, floorcap, dfloor_flux_limit, FOFC, update_mcap, merge_lat, wellbalance, the sponge and dt_retry are **not** in the ledger: none of them books into it.

## Binary
- **LZC** = branch **`caltech-lzclass-920f2644d`** (pushed) = vet-conserv-lzab2 b5700b17c + 920f2644d. The new key is `<hydro>/p_ledger_lz` (read only when named; 3-D spherical polar with `<rad_vet>`). It turns on p_ledger with L_z:
  - **Gas + radiation L_z**: `Lz_tot` less the ledger, radiation from the x3 face F. It is booked by the p_ledger classes `hflx ntel src rad cut fwb flr mrg gap`. `hflx` is the gas x1 torque (am_bflx); cut, fwb and fref stay in `rad`.
  - **`lzb`**: the ledger's growth.
  - **`rbt`**: the part of `lzb` that is the radiation's torque through the physical x1 faces, booked like LzBook.
  - Output: `CLEAN LZLEDGER` lines; `user.hst` gets `LzL_tot LzL_hflx ... LzL_gap LzL_lzb LzL_rbt` in place of `Pr_*`.
  - It only reads the state.
- Orion CPU gate, job 243739:
  - star_vet md5 ac6824fc, tiny_bsg_rel5 pass1 3.720e-11 PASS.
  - Suites: quick 28/28; targeted 7/7 (conservation_lz, floorcap_spin, conservation_sp/cs, restart_sp/cs, p_ledger).
  - p_ledger test: lz_end 3.4e-16, lz_bitwise 0, lz_restart_lines 0.
  - Tiny bsg 20 cycles, on vs off: dump data, hydro.hst and the 32 old user.hst columns are identical. Restarted at 10, user.hst including LzL_* is bitwise and LZLEDGER lines differ 0/10. 0 CLEAN WARNING.
- Build as REL: CUDA 12.9 + g++13, star_vet, never CUDA 12.2. Tiny gate first (pass1 3.720e-11); post md5 and prov.

## Orion CPU class table (tiny bsg, 20 cycles, rate x dt summed over the run, g cm^2/s)
```
hflx -1.19e33  ntel 4.83e42  src -4.83e42  rad 4.125e34  cut 0  fwb 0  flr 9.3e23  mrg 6.18e32  gap 0
lzb -4.1245e34 (= CLEAN LZB growth)   rbt -4.1163e34 (99.8 % of lzb)   END sum 1.6e26
```
- ntel and src cancel, as the flux divergence and the coordinate/am sources should. rad + lzb = 7e30. mrg (merge_lat average, 6e32) and hflx are the gas's own L_z changes outside the ledger.

## Arms (4 single-GPU arms on one H200 node; fresh from the LZC templates `inputs/stars/star_bsg.athinput`, `star_t9l62.athinput`; 100 cycles; `rad_vet/clean_verbose=1`, `output1/dcycle=1`)
| arm | extra keys | question |
|---|---|---|
| LZC-bsg | `hydro/p_ledger_lz=true` | class table on the production bsg |
| LZC-t9 | `hydro/p_ledger_lz=true` | the same on t9l62 (27zf: 3.5e38 per 300 cycles) |
| LZC-m35 | `hydro/p_ledger_lz=true` (star_m35mms) | third star (27zf: -2.7e34) |
| LZN-bsg | none | byte identity with LZC-bsg: hydro.hst, the old user.hst columns, final hydro_w data; and the cost (s/cycle) |

## Report (NOTE caltech-orion-lzclass with this letter; files to caltech_results/LZCLASS/<arm>/)
1. Per arm, per cycle: the LzL_* columns and the CLEAN LZLEDGER line. Sums over the arm (rate x dt) per class, plus lzb, rbt, rbt/lzb, rad + lzb, and END max/RMS.
2. CLEAN LZB ledger first/last (growth) against the sum of lzb.
3. CLEAN WARNING / FATAL / DT_RETRY counts; s/cycle (LZC vs LZN, same node).
4. LZN vs LZC byte identity (as above).

## What we expect and the fix we will propose
- Expect rbt/lzb ~ 1 on every star. The ledger is then the radiation's L_z leaving through the outer face: a real discrete boundary flux, not a lost-L_z error. The radiation's counterpart of Lz_hflx is missing from the budget.
- Proposal, accounting only, no result change: book rbt as a radiation boundary flux (an `Lz_rflx` history column next to Lz_hflx), so that dLz_tot/dt = -(Lz_hflx + Lz_rflx). The ledger then holds only lzb - rbt (0.2 % on tiny bsg).
- Physics (a separate decision; it changes results): D_{r phi} != 0 at the open outer face of a non-rotating star comes from the closure's angle set (the one-sided `pe` at the open/Marshak face takes the outermost cell's D_{r phi} E). A face closure at the open x1 face consistent with the boundary is one option: D_{r phi} = F_r F_phi / (c^2 E) there, which is 0 for F_phi = 0. A phi-symmetrised angle set is another.

-- CONSERV-LZCLASS (orion)
