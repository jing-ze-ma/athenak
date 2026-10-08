# NOTE 2026-10-08 Caltech: He giant N445 done, N897 remap + smoke PASS

From: Caltech. To: viper, DeltaAI. Docs only.
- N445 1-H200 run (4212784) reached tlim 9.004208e5 (10.42 d) in 13.1 h, 0 FATAL; final rst hegiant.00021.rst (t 900420.8, cycle 47862).
- Remap (he_remap_giant.py, recipe 1307089d) -> N897 hegiant.00021.rst (2.27 GB): rad_cons_* 4-7e-16, rad_Fr_change_rel 0.72 %,
  rad_hyd_bad_after 0, rad_F_ratio_max 1 (= viper's expectations).
- N897 needs >= 2 H200: on 1 H200 Kokkos failed to allocate 119 GiB.
- Smoke 4238905 (2 H200, 2 ranks, viper first-link keys, 60 cycles): rc 0, 0 FATAL/NaN/NON-CONVERGED; dt 0.024 (first step after
  refill) -> 18.58 (viper 9.22 -> 18.59); Picard mean 16.9 max 30 (viper 18.5 / 42); 1.5 s/cycle on 2 H200
  (~38 h on 2 H200 from 10.42 d to 30 d at dt 18.6).
- N897 production NOT started on Caltech: the user decides (DeltaAI link 1 3337674 and viper hegiant897 are the other copies).
