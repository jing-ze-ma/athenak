# m1-fast5-box RESUME (2026-09-25 night)

Branch m1-fast5-box (from rt-integration 2c2a0d78), worktree /viper/ptmp2/jinma/wt_fast5box.
Run tree /viper/ptmp2/jinma/fast5box_0925 (bin/, runs/, scripts/, smoke/, off/, gates_*).

## State
- Lever 1 DONE (WIP commit 2271b0ce + this commit): `<rad_m1>/implicit_precond = mg_gf`
  (keys implicit_gf_modes = 2, implicit_mg_levels = 3), opt-in; default and mg untouched.
  Code: rad_m1_precond.cpp (ImplicitGFInit/Build/Pre/Add, GFLayerReduce/Sum/Project),
  rad_m1.hpp (gf_* members), rad_m1_implicit.cpp (+choice mg_gf, key parsing, ~12 lines).
  Method: 25 lowest (x2,x3) Fourier modes per x1 layer; coarse op = mode-diagonal Galerkin
  of the layer-mean stencil (pentadiagonal in x1), LU + inverse on device once per pass;
  applied coarse-first in projection form z = P Ac^-1 P^T r + M_mg (r - Pi r).
  Offline (off/, study3.py on coarse2 box dump): mg3 9 it -> 5 it.
- GPU box timing (runs/t5, job 11980650, binary n5, 3 reps, cycles 10-50):
  1 GPU mg3 44.11 -> mg_gf 41.65 ms/cycle (-5.6 %), it/solve 8.5 -> 4.5;
  2 GPUs 50.15 -> 49.28 (-1.7 %), 8.7 -> 5.4. Timers (runs/p5/tm_*): Krylov 13.70 -> 11.26
  ms/cycle, 0.795 -> 1.19 ms/it, radiation/hydro 2.08 -> 1.95. gf23 is best of gf13/gf22.
- n6 = n5 + 4-accumulator sum kernel (committed source = n6).

## Gates
- PASSED (binary n2, gates_eval.txt): default + mg3 CPU bitwise base vs new (12/12);
  mg_gf vs rbgs_fwd F1top 1.4-1.6e-10 box / 3.1e-9 slab = level of mg3 vs rbgs_fwd in the
  same run (1.5e-10); mg_gf vs mg3 vet_sc PASS; 1 vs 2 ranks PASS; 0 non-converged.
- PASSED GPU: mg3 rst + hst bitwise vs base at 1 and 2 GPUs (runs/t1/b?base_1 vs t5/b?mg3_1).
- PENDING (running at handover): scripts/cpu_gates.sh n6 -> gates_eval_n6.txt,
  rstgf_n6.log (restart bitwise, mg_gf 2 ranks, Eddington + vet_sc), tst_rad_m1_n6.log;
  GPU timing job 11980862 -> runs/t6 (python3 scripts/tsum.py runs/t6 10 50), and check
  cmp runs/t1/b1base_1/rst/m1slab.00001.rst runs/t6/b1mg3_1/rst/... (same for b2).
- Not needed: radwave order tests (default path bitwise; mg_gf is a preconditioner only).

## Next steps
1. Read the pending gate files; if all pass, squash-describe the commit (no merge/push).
2. Lever 2 (not started): ~1.3 ms/cycle GPU idle from the two host syncs per BiCGStab
   iteration (idle before m1_gf_red / m1_impl_bcg2_upd, runs/p5/prof/*.lsum);
   implicit_krylov_dev refuses mg/mg_gf and is 1-rank only - making M1PCRXD/GFProject read
   beta/omega from the device kd array is the path.
3. mg_gf remaining cost: per pass sbar+fac+sum ~0.3 ms (could build once per solve),
   per application red 26 + solve 24 + rp 13 + pro0 +9 + h 5 us; 2 GPUs pay one host
   MPI_Allreduce per application.
4. Lever 3 (SC sweep 4.1 ms, stencil build m1_impl_stb 0.52 ms/call with c[19] scratch)
   not started; fast4 already tried the stb unrolling.
5. Delete src_n1..n6 and b_* dirs in the run tree when done (build dirs already removed).
