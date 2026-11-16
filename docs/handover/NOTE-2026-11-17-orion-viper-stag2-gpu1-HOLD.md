# NOTE 2026-11-17 Orion -> viper (USER): HOLD the stag2 GPU test (TASK e38e2b0de) until the follow-up package

USER decision: please do not start TASK-2026-11-17-orion-viper-stag2-gpu1 (branch clean-cc-stag2-gpu1 9689f4216) yet. It went out early,
before these items were in:
1. **CUDA compile:** rad_m1_clean.cpp, which holds most of the new staggered kernels, has not finished compiling under nvcc on Orion
   (cicc > 45 min). Its GPU readiness is unverified.
2. **flz q wall regression:** in the integrated set the wall gas velocity points outward again (+0.003..+0.010; cell -0.0032).
3. **F-solve with the transverse coupling in its matrix** (clean_stag_fmat = mg): F is converged to 1e-12 with the scores unchanged. It is
   green on Orion and not yet in the branch.
4. **Cubed sphere:** the transverse kinetic stress is being dropped on cs, to fix rhse_cs WB (2.96e-9 against 9.8e-9 now, 1.5x faster).
5. **Marshak front fixes:** the beta-term consistency fix 7c9bc150e.

Orion will send a follow-up TASK on a new branch hash with the same decks, refreshed CPU references and gates. Nothing in the cell layout
or in viper's code is affected; the branch is bitwise to 0678599eb with staggered off.
