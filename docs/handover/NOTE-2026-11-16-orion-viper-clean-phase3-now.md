# NOTE 2026-11-16 Orion -> viper (USER): clean scheme handed over NOW (phase 3 immediately)

USER: "hand it over to viper now". From this NOTE on **viper owns origin/clean-cc**: code, merges, pushes, CPU and
GPU testing. Orion pushes nothing more to origin/clean-cc and only helps on CPU when viper or the user asks
(QUESTION/TASK files here). The rest of HANDOVER-2026-11-16-clean-scheme-orion-to-viper.md (scheme, keys, binding
rules, ranked open problems) stands; read its sections 1 and 6 as superseded by this NOTE.

## State
- origin/clean-cc = **62d71accc** (push 13 = push 12 WB limiter fix + h-free Rhie-Chow / time order 2; NOTE ed0ec221c).
- Your fix 073378e14 (clean-cc-gpufix-viper) is NOT on origin/clean-cc. Orion merged it locally onto 62d71accc
  (clean merge, no conflicts) and ran the full CPU regression: **bitwise to push 13 on every deck** (25 decks,
  every bin/tab/hst), restarts bitwise, zero CLEAN WARNING. Merge it and push when your GPU G3/G4 reruns pass.
- G1 cb score flag: CPU run.log implicit_beam lines are identical to your GPU ones (ang = -45 / 45, same rays and
  patches), so the 0.4707 vs 0.4600 difference is NOT the log angle; with identical dumps it has to be the scorer
  environment (Python/numpy version or the dump file read) — check on your side.

## Orion work handed to viper (take, adapt or drop)
1. Clump / porous-photosphere test pgen (`m1_test = clump`, 6 files, +425): handover_1116_clean/
   clean-cc-clump-2db207987.patch (`git am` onto 62d71accc). Results (SC best, Ftop +2.5..+8.8 %; FV arms need a
   rerun on head) are open problem 8 of the handover.
2. Regression harness: handover_1116_clean/reg/ (follow-up commit from Orion: decks, cmds, analysis, restart harness,
   WB check, reference numbers for 62d71accc, README). Run it before every push.
3. Member-initialisation audit (Orion worker, in progress): every RadiationM1 / clean member set only on some paths
   but read on others, and 2-D last-face out-of-range reads of the 073378e14 kind. Result comes as a patch in
   handover_1116_clean/ with a NOTE; your A100 rhse_sp2 eint <= 0 may be in this class.
4. Open user decision: time-dependent rays (open problem 4; Orion recommended implementing for FV/lld first, no
   rotation). Ask the user.
5. IDORT suite page (artifact WTbarCvoRbP4oNhmTGCA4f, built on 3609448c4): its time-order maps predate push 13;
   rerun is yours if wanted (pipeline: Orion /orion/ptmp/jinma/clean_1115/idort_page, ask for a copy).
