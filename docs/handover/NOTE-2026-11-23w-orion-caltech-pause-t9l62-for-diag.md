# NOTE 2026-11-23w orion -> Caltech: USER asks to PAUSE IA11 t9l62 and run the 23q diagnostic on its node now

USER (11-23 ~14:45 PST): "ask caltech to stop t9 for a while to do the diagnostic".

## 1. Pause IA11 t9l62 cleanly
- `touch prod_ia11_1123/t9l62/STOP`. Let it finish its cycle and write the final rst, and confirm "Terminating on stop request".
- Keep that rst (and all of t9l62's hst and logs). Post its name, t, cycle and md5.
- Leave the bsg and m35mms IA11 runs and holder 7333072 alone.

## 2. Run 23q arm A on that node (hpc-sm-02-09, 4 H200)
- Binary: the IA9 build `athena_ia9_77551e823` (md5 bb0f2fbe). The rst is IA9's, so do NOT use IA11.
- Restart from `prod_ia9_1123/t9l62/rst/t9l62.00001.rst` (t = 4.757564e5) in a NEW run dir, with the same layout as the IA9 run (4 ranks, 2 blocks/GPU).
- On the command line add only these diagnostic keys (they do not change the state):
  `hydro/dt_report=1 hydro/floor_audit=1 time/nlim=300 output9/file_type=log output9/dcycle=1`
  - output9 is a new block number, so it does not collide with the rst's stored outputs.
  - hst is already every cycle.
  - No rst / bin outputs: set their dt to 1e30 in a deck copy if needed.
- Arm B (xpat) is cancelled (23t). Only arm A.
- Expected: ~300 cycles in a few minutes. dt is small in this phase, so it may take ~10-20 min.

## 3. Post, then resume
- Post log.0(.gz) + user.hst + hydro.hst + run.log into caltech_results/IA9_T9DIAG/A/ (gzip). No rst.
- Then RESUME IA11 t9l62 from the rst written in step 1: same binary athena_ia11_bd9553537, same deck, same launcher, and nothing changed on the command line except -t.
- Post the first cycles of the resumed run so we can check bitwise continuation (hst rows after the restart).
- Delete the diagnostic run dir's dumps after posting.

Why: this is the only way to see, on the exact full-resolution state, where the ~3e42 erg went in IA9's floor_pay cycles. It decides our next binary (IA12), which we will send separately.

-- orion (coordinator)
