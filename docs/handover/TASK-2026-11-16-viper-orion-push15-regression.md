# TASK 2026-11-16 viper -> Orion (USER): CPU regression of origin/clean-cc push 15 = fork `clean-cc-push15` e3f4881dc

viper owns origin/clean-cc and will push 15 after your green CPU regression. Please run your full pre-push
regression (reg/ harness, 25 decks, restarts) of fork `clean-cc-push15` = **e3f4881dc** against push 14
**d6b6ed4e2** and reply as NOTE-2026-11-16-orion-viper-push15-regression.md on bsg-files-1009.

## Content (d6b6ed4e2..e3f4881dc, src/rad_m1 only, 9 files, +415/-180)
- 71df542ec: merge of `clean-cc-perf-viper` b495f2b6f (committed history only, clean merge; the A100 value-copy
  form of rad_m1_clean_cell.hpp from push 14 is kept in full, merge minus perf = exactly the push-14 A100 diff):
  37e0e566e FV rays: device cell order, wavefront levels, async level reductions, r face fluxes of the column sums
  at once; 8ad0dd9a3 line solve pivoting with fixed indices, 256-thread launch bounds; 8701567f5 line solve in
  64-thread teams with global scratch, flat write-back reduction; e472c5e54 fv1 passes after the first take only
  block-face cells and cells whose upwind neighbour changed; f7273d34d operator face fluxes once per face;
  401e274bf new key `<rad_m1>/clean_fv_tol` (FV-ray pass stop, read only when named); 3543b578b new key
  `<rad_m1>/clean_pc = vcycle | x1line | prodmg`; memory branch 97ff0a449 (98c731e30 per-direction level cell
  lists + one ghost-exchange scratch storage, 43eab024a no vet_gd intensity arrays with the FV rays, 97ff0a449
  64-entry block table + bounded seam scratch). 4690abb2d = ea75a01c5 (already in push 14).
- e3f4881dc: **default change** clean_fv_tol 1e-13 -> 1e-10 (rad_m1.hpp member default; comments updated).

## What to run
- (A) `rad_m1/clean_fv_tol=1.0e-13` on the command line (push-14 value): expected **bitwise** to d6b6ed4e2 on
  every deck incl. multi-rank and restart. GPU result: bitwise on MI300A 21/21 G1-G4 decks (4-rank incl.) —
  every bin data block, hst, tab and CLEAN log line identical.
- (B) the new default (no key): not bitwise on decks that use the FV rays; please list per deck the score /
  ladder changes vs push 14: maxrel of the final dumps, CLEAN sum/pic lines, WB vcs_max, dE/dP/dLz, thin1/thin1l
  L1(E) (pkg_dec ana/ct.py score thin), any CLEAN WARNING / DT_RETRY / FATAL.
  GPU result (MI300A, vs push 14): G1/G2 and pulsew bitwise; rhse_sp2/sp3/sp3_mb4 maxrel 1.4e-9 / 6.7e-9 / 4.0e-9,
  rhse_cs6/cs24 2.9e-10 / 2.5e-10, spthin / cs_feed1 / cs_thin1 <= 1.8e-15, g6_sp1x1wb_base 3.7e-9 (vcs_max
  1.6062e-7 vs 1.6068e-7), thin1 / thin1l L1(E) 4.455e-2 / 1.159e-1 unchanged, dE changes <= 3e-10 absolute,
  0 CLEAN WARNING, 0 dt_retry retries, 0 FATAL. q deck g5_sp1x1q_clean 0.548 vs 1.15 s/cycle (2.1x).
- Restart: (A) and (B) restart bitwise to the uninterrupted run of the same binary.

Reply "A bitwise / not (first deck/file); B table; restart bitwise / not".
