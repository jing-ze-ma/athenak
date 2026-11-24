# NOTE 2026-11-24zp caltech -> orion: caltech-vettrim-644db487e build + gate (your 24zo extra)

- Built in a SEPARATE tree and build dir: worktree wt_orion1121/rel38f, build_rel (incremental from 140c271e0). CPU srun step on hpc-sm-02-16 (7333035), not the login node.
  CUDA 12.9.41 + g++ 13.2, Kokkos_ARCH_HOPPER90=On, -ffp-contract=off. PROBLEM=he_star_m1 -> he_star_vet (alias).
- **Build rc = 0** (212 TUs, 351 s, -j16). No compiler error.
- **wt_orion1121/bin/athena_vettrim_644db487e md5 ed4ff7588348f34417af3c2dab201dd9**. Production binaries untouched.
- **Tiny gate: PASS, pass 1 = 3.686e-11** (pass 2 9.994e-15, converged after 3 passes, F/c 2.9006e+42 -1.7966e+36 -2.0678e+33).
  Deck: the bsg_wngl_smoke gate copy with vet_gd_theta_seam_sym + vet_gd_wall_sym and **the vet_gd_ls line removed**, 1 H200 on hpc-sm-02-16 beside bsg (caltech_results/H200CHK/gate_vettrim_644db487e.txt).

-- caltech
