# TASK 2026-11-16 viper -> Orion (USER): CPU regression of fork `clean-cc-push16` 51e2ce016 vs push 15 dfe21a746

viper owns origin/clean-cc and will push 16 only after a green CPU regression from you (plus the MI300A and A100
gates). Please run your full pre-push regression (reg/ harness, all decks, restarts) of fork `clean-cc-push16` =
**51e2ce016** against push 15 **dfe21a746** and reply as NOTE-2026-11-16-orion-viper-push16-regression.md on
bsg-files-1009.

## Content (dfe21a746..51e2ce016, src/rad_m1 only, 9 files, +1193/-110)
- 57a08c09f (= f0fd11874, bug fix, always on): the clean scheme's gas internal energy is the total less
  d phi (etotgrav) and |B|^2/2 (MHD) through CleanGasEint.
- 0a3df762d (= 599b781b8, bug fix, always on): Newton change of E over max(|E|, |E_old|).
- ff74a049e (= d7362e956): keys `clean_gas_exact`, `clean_yres`, `clean_nl_norm = global | cell`,
  `clean_nl_efloor` (added opt-in, defaults below).
- d589df2d7 (r136cost): `clean_pc_factor` (x1line keeps its LU factors, bitwise), `clean_pc_reuse`,
  `clean_conv_est`, `clean_nl_print`; the x1line preconditioner no longer builds the coarse MG levels.
  Merged with dfe21a746: both x1line kernels use teams of one thread on the host.
- b7964ecaa (rot2): `clean_rot_ring = off | live | jump` (rotating vet_gd FV union as a ring; restart channels).
- 44c0fedce (fv1pass): `clean_fv_joint` (fv1 levels over all MeshBlocks of a rank, one exact pass; falls back
  to the per-block levels with sp, multilevel, rotation or any cell reading another rank) + pass-flag realloc.
- 51e2ce016 **defaults**: clean_gas_exact = true, clean_yres = true, clean_pc_factor = true, clean_fv_joint =
  true, clean_rot_ring = live when not named and the vet_gd FV rays rotate in time with vet_rotate_advance,
  _align and vet_gd_rotate_xfade. Opt-in: clean_pc_reuse, clean_conv_est, clean_nl_norm = cell.

## What to run
- (A) push 16 with the new defaults forced off on the command line:
  `rad_m1/clean_gas_exact=false rad_m1/clean_yres=false rad_m1/clean_pc_factor=false rad_m1/clean_fv_joint=false
  rad_m1/clean_rot_ring=off`. Expected **bitwise** to dfe21a746 on every deck incl. multi-rank and restart,
  except decks that have etotgrav or MHD with the clean scheme (the f0fd11874 fix); list those separately.
- (B) the push-16 defaults (no keys): per deck vs push 15: maxrel of the final dumps, CLEAN sum/pic lines,
  WB vcs_max, dE/dP/dLz, thin1/thin1l L1(E) (pkg_dec ana/ct.py score thin), any CLEAN WARNING / DT_RETRY /
  FATAL. Please include one x1line deck (`rad_m1/clean_pc=x1line`) in both arms.
- Restart: (A) and (B) restart bitwise to the uninterrupted run of the same binary. Note the ring: a rotating
  deck restarted from a push-15 file is expected to skip the slope record (channel count differs: first step backward Euler, WARNING).
- UPDATE: test the branch head **a97e368cb** (51e2ce016 + c249955c4 nvcc-only forms of two
  lambdas + a97e368cb = race fix dad6ddabd from clean-cc-race-viper: each x1-line team owns its
  rows in cln_lw instead of team scratch; host teams of one thread kept). Expected bitwise to
  51e2ce016 on CPU.
- UPDATE 2: test the head **0aee9a208**: clean_yres is back to opt-in (default false). With it on,
  the first backward-Euler step of rhse_cs6 / rhse_cs24 does not converge (CLEAN WARNING +
  dt_retry, MI300A and A100; bisect: off -> clean, gas_exact or fv_joint off -> still warns).
  So for (A) add nothing new (rad_m1/clean_yres=false is now the default); (B) = no keys.
- UPDATE 3: final candidate head **0678599eb** = 0aee9a208 + your init audit (4cacdae21,
  8db6f6141) + steady-start predict fix (5e33f056a) + your dt_retry 0001/0002 (d320347b3,
  4ac0c0e49) + clean_fv_joint back to opt-in (GPU: cs_thin1 +38 % A100 / +62 % MI300A).
  Defaults on: clean_gas_exact, clean_pc_factor, clean_rot_ring = live (rotating sets).
  (A) keys: rad_m1/clean_gas_exact=false rad_m1/clean_pc_factor=false rad_m1/clean_rot_ring=off.
- UPDATE 4: final candidate head **dd468ab99** = 0678599eb + 0a03129a8 (clean_rot_ring default in
  the constructor: the default set in CleanInit broke restarts of rotating runs) + c40170eab (your
  cubed-sphere FV-ray default in the constructor) + 665ae99de (FV-ray default before the ring
  default) + dd468ab99 (comment within 90 columns). Expected bitwise to 0678599eb on fresh runs.
- GPU (MI300A apudev, A100 Raven) gates are running now; results follow in a NOTE.

Reply "A bitwise / not (first deck/file); B table; restart bitwise / not".
