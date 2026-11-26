# NOTE 2026-11-27zl orion -> Caltech: scpr fixes (your 27zh reply), PR timing rerun + accuracy script

Thanks for 27zh. Two errors were mine:
- The PR fatals on a command-line key that the deck does not name: `pr_sph.athinput` now has
  `iter_min = 1` (pull scpr-1127 again).
- `sph_err.py` assumed reader keys that do not exist (`mb_x1v`), and the PR names its J `rad_J`.
  Both are fixed (cell centres from x1min/x1max, `sc_rad_J` maps to `rad_J`).

Please, on the same node/binaries (PR md5 a77bc7ee, ours rel6 suite 8cc98fb3):
1. Rerun **G1_P1, G1_P2, G1_P1, G1_P2 and G4_P1** of 27zh unchanged (with the timer). Post the s/cycle,
   `cpu time used`, the `Iteration stops` lines and the top 25 lines of kp.0.
2. Run `sph_err.py` on the A64/A32 bins if you still have them. Otherwise rerun A64_pr, A64_our,
   A32_pr, A32_our (a few seconds each) and then the script: `sph_err.py <bin> sc_rad_J` and
   `sph_err.py <bin> vet_vj gl 8 20`. Post the printed lines (interior error + exterior profile).
3. If 27zk (thick sphere chi = 50: 5a ours Cartesian, 5b ours sp with `our_sph_sp.athinput`, 5c PR)
   has not run yet, run it now with the fixed script and deck.
Then delete the bins. Reply caltech -> orion "scpr3".

-- orion (SCPR, for the coordinator)
