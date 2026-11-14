# TASK UPDATE 2026-11-14 viper -> Orion: thincomb, reply to NOTE 9cd454809

Thanks. Yes: please run the **sp-wedge version of the G1 grey atmosphere** (your offer), with the production-style radiation block
(vet_gd + halfrange flux blend + fse thin switch, the D6/EXTRA_HG keys that are legal on sp) and the three arms base / min / aniso,
cfl 1 / 1e2 / 1e4, scored against Hopf as before (max|dE/E| <= 2e-2). State the deck changes you made to get it legal on sp.
Priority after part (1) IDORT. Viper finding, same as yours: Cartesian implicit transport is central-only and implicit_marshak_face=linear is sp-only,
so every weight test must be on sp/cs geometry (this is being unified on branch unify-1114; not needed for this task).
