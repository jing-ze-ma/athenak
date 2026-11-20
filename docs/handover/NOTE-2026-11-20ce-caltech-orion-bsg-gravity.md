# NOTE 2026-11-20ce Caltech -> orion (USER): BSG gravity in the earlier Linhao / Ma+2026 reproduction (answer to 77b95b844)

Checked from the decks on Caltech (bsg_1001, bsg_1009, bsgprod_1012, bsgli_1014, bsgli_1104, bsgma_1104, bsgsp_1103).

1. **GM:** every reproduction deck used **he_gm = 2.6542488e27 = G x 20.00 Msun**, **not** his 2.62502e27 (19.78 Msun).
   This covers the true repro bsg_1009/bsg3d_truerepro2, the 39-d run bsgprod_1012/bsgP_seed_fin3, the Linhao-exact family
   bsgli_1104/bsgLI_*, and bsgma_1104. The bsgli_1104 deck header says why: "he_gravity pointmass, he_gm = G 20 Msun |
   [PAPER] (his GM value unknown)". We took the paper's M = 20 Msun because we did not have his GM at the time.
   The one exception is bsgsp_1103, our own BSG setup and not a reproduction: he_gm = G m(20 Rsun) = 2.6309351439e27 (19.82 Msun).
   That came from our MESA-based column (tools/ic_star mode mesa), with he_gm_column = true.
2. **mu:** **0.62** (the paper value) in every reproduction deck, with an ideal gas, gamma 5/3 and kappa_es const 0.34. 0.614 was never used.
   bsgsp_1103 used our own EOS, not ideal 0.62.
3. **Gravity form:** a **point mass only**. There was no envelope self-gravity and no g(r) from the profile.
   - The reproduction decks have he_gravity = pointmass, except bsgma_1104, which has he_gravity = wb (well-balanced, wb_option polytropic)
     around the same G x 20 Msun point mass.
   - The only m(r)-type gravity was bsgsp_1103: he_gm_column = true, so G m(r) = G m(r_bot) + G int 4 pi r^2 rho dr of our column.
   - The IC atmosphere above the profile (r > 52.49 Rsun) used he_ic_linhao_gext = invsq, i.e. HSE in GM/r^2. The prod_hisg variant used r2,
     which is his g ~ r^+2.
4. **HSE check of the Linhao profile against GM:** we have no record of one in the reproduction work. The profile was taken as given at
   mu 0.62 / 20 Msun. Our only HSE audit was for bsgsp_1103, our own column (SETUP.md: force audit; deep HSE drift < 1.3e-4 below 0.5 R_ph).
   That is a different IC, so it says nothing about 19.78 vs 19.9.
   Your numbers fit this: at mu 0.62 the profile's HSE mass, 19.68-19.94 Msun, brackets both 19.78 and 20.0. Our 20.0 sits near the top
   of that range, so our point mass was 0-1.6 % strong over the domain where his 19.78 was 0.2-1.1 % weak at mu 0.614.
5. **Where it is recorded:** the deck header table in bsgli_1104/bsgLI_exact.athinput, plus bsgprod_1012/bsgP_seed_fin3.athinput
   ("he_gm = 2.6542488e27  # G M, M = 20 Msun (point mass)"). For bsgsp: bsgsp_1103/SETUP.md and ic_star_blend/README.md.
   There was no separate NOTE; it was never discussed as a decision.

-- Caltech
