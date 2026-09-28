---
name: correlated-k-shared-module
description: The Exo-FMS correlated-k opacity layer extracted from the hot Jupiter pgen into src/utils/correlated_k.hpp (2026-09-07), so red_giant and others can use it; includes the Rosseland builder that couples it to radiative conduction. Gated bitwise
metadata:
  type: project
---

**Done 2026-09-07 on orion at the user's request** ("extract the ck rt into a shared
module" ... "make sure it supports merging with conduction"). UNCOMMITTED.

## What moved
`src/utils/correlated_k.hpp`, namespace `correlated_k`, 1071 lines. From
`src/pgen/deep_hot_jupiter_rt.cpp` (6082 -> 5160 lines): the OPACITY layer (old lines
536-1408: k-table reader, CIA/Rayleigh/H- continuum, FastChem chemistry table,
band-integrated Planck fractions, `ck_tp_index`/`ck_kappa`/`ck_continuum`/
`ck_planck_bands`, `ck_selftest`, `ck_rt_selftest`) AND the Rosseland builder (old
3140-3238).

**The two-stream SOLVER did NOT move.** It is written around the hot Jupiter's
irradiation and arrays. Stage 2 if a second pgen needs band RT rather than just band
opacities.

## Merging with conduction (the point of the exercise)
`ck_build_rosseland_table(Conduction *pc)` (was `(Mesh*)`) tabulates the Rosseland mean
of the SAME k-table + continuum on the table's own (T,p) grid and installs it as
`pc->rad_kr_tab` etc. That is what makes the tau blend consistent: two-stream above,
diffusion below, ONE opacity, so the handover does not jump. Its old `if (!rt_ck) fatal`
became `if (ck_lk_ptr == nullptr) fatal`, which is the module-appropriate test.

## Design decisions
- The module owns NO solver state. Three signatures now take it as arguments:
  `read_ck_table(file, pcut_bar)`, `build_planck_fractions(pcut_bar)`,
  `read_ck_continuum(dir, swfile, star_teff)`. Previously they read the pgen globals
  `rt_ck_pcut`/`rt_star_teff` (pcut only for an advisory message; star_teff genuinely
  builds the stellar band fractions).
- The pgen keeps every use spelled as before via a block of `using correlated_k::X;`.
  Total ADDED to the pgen: the include, a 3-line comment, ~45 usings, and 6 lines at the
  four call sites. Nothing else.
- `inline` on all namespace-scope definitions. TRAP: `KOKKOS_INLINE_FUNCTION` already
  expands to `inline`, so prefixing the following line gives "duplicate 'inline'" on the
  5 device functions.

## THE GATE (do this again for any further move)
Baseline `problem/ck_dump_file` column from the pre-refactor binary, then compare:
`/orion/ptmp/jinma/Athenak/ck_extract/{base,after2}`. Result: column BITWISE IDENTICAL
(md5 ac3aac695ec76c7cbc1a330a889f694c), history identical, self-test and Rosseland
startup reports identical. `test_rad_dhj_ck_cpu` PASSES (44 s).

## Style
The header is cpplint-clean at 90 columns. It inherited 22 over-long lines from the pgen
(which itself has 139 and 161 trailing-whitespace lines, both PRE-EXISTING at HEAD and
untouched); those were re-wrapped in the new file only. Two fatal messages that named
`deep_hot_jupiter_rt` now name `correlated_k`.

## Next
A PLANCK mean built the same way. The Rosseland mean is a poor stand-in for it in an
optically thin relaxation, which is exactly what `red_giant` currently does -- see
[[red-giant-envelope-project]]. The weighting is the only difference.

## STAGE 2 (same evening): the two-stream SOLVER is out too
`src/utils/two_stream_rt.hpp` (1624 lines) and `src/utils/atm_column.hpp` (105). The pgen
is now **3648 lines, from 6082** before any of this.

- `two_stream_rt`: the RT macros (RT_NB/RT_NNC/RT_CACHE/RT_FP32, RtF), all 24 rt_*
  config + scratch globals, `par_reduce_clip3/4`, `LimitRTSource`,
  `RTSourceLimiterWarn`, `get_albedo`, `get_picket_fence_coeff`, `get_kapr`, `get_Tint`,
  and `picket_fence_two_stream_RT` (1108 lines: the ck AND grey picket-fence sweeps).
- `atm_column`: the small pure helpers its kernels call -- `CSCellAngles`, `GravAccAt`,
  `GravPotAt`, `TideAccR/T/P`, `EffGravAt`, `TGuess`. Put in their OWN header because a
  gravity model does not belong in an RT module.
- Configuration still flows through the namespace-scope rt_* variables and
  `pm->pgen->hot_jupiter_param` (a generic carrier despite the name: Teq, grav, ap,
  omega, Rgas, met + flags; a self-luminous object sets Teq = 0). Turning that into an
  explicit config struct is the next tidy-up and was deliberately NOT done in the same
  change, so the move could be proved byte-for-byte.

### Four traps, all found by the compiler
1. `body_end` brace counting must STRIP comments and string literals first, or `{` in a
   comment desynchronises it (this failed on the first function).
2. The pgen's `namespace {` opener sat INSIDE the moved global block, so it went to the
   header and left an orphan `}  // namespace`. The anonymous namespace wrapped ONLY
   those globals, so the fix was to delete both, not to re-open it -- re-opening it at
   the top put forward declarations and definitions in different scopes and every call
   became ambiguous.
3. Deleting a forward declaration leaves its `KOKKOS_INLINE_FUNCTION` line behind; five
   of them stacked up and gave "duplicate 'inline'" pointing INSIDE Kokkos_Macros.hpp.
4. The comment-aware unit walk gave `LimitRTSource`'s 27-line doc to `CSCellAngles`,
   because in that file the doc and the body are separated by another function.

### The gate, again
Same baseline column. After BOTH stages plus every style fix: column BITWISE IDENTICAL,
history identical, `test_rad_dhj_ck_cpu` passes (40 s). Both new headers are cpplint-clean
at 90 columns with no trailing whitespace; the pgen's long-line count fell 162 -> 109 and
its 161 trailing-whitespace lines are pre-existing and untouched.

### Still to do
red_giant does not yet CALL any of this -- that is the point of the exercise and the next
step: include two_stream_rt.hpp, fill hot_jupiter_param with Teq = 0, allocate the rt_*
scratch, and replace the grey Eddington relaxation with the band solver.
