# tests_bench_audit_0922 -- NOTES.md / README.md vs. actual .athinput audit (2026-09-22)

Read-only audit. Nothing under `bench/` was modified.

## Method

1. Enumerated every `NOTES.md` / `README.md` / `notes.md` directly inside a `bench/*`
   directory (depth <= 2): **58 files**.
2. Built the set of every parameter name AthenaK's `ParameterInput::Get*`/`GetOrAdd*`
   actually reads at this HEAD (`grep -rhoE` over `src/`, 992 distinct `block/name` pairs).
3. Scanned each notes file for lines of the form `` `name` = value `` / `name = value` /
   `name on|off|true|false` where `name` is one of those known parameter names (>=5 chars,
   contains `_`, to cut noise from physics symbols like `a`, `f`, `rho`, `GM` that also
   happen to be legal keys for other pgens) -- 41 candidate claims after that filter.
4. Compared each claim against the `.athinput` in the exact sub-directory the claim's own
   context names (arm directory, e.g. `nowb/`, `etalow/`, `mhd_thuni/`); fell back to the
   shallowest `.athinput` in the notes file's own directory only when no arm was named.
5. Manually resolved every automated MISMATCH/NOKEY hit by reading the surrounding
   paragraph and the named arm's actual input, since the automated first pass produces
   false positives whenever a comparison table names an arm the naive matcher didn't
   pick (most of the 41 raw hits were exactly this: real per-arm overrides, correctly
   set in the arm's own file, that only *looked* wrong against a sibling arm's file).
6. Cross-checked every explicit `md5 <hash>` claim in a `BUILD_COMMIT.txt` / `NOTES.md`
   / `README.md` against `md5sum` of the actual binary, where the binary still exists.

Coverage: all 58 notes/readme files were scanned by the automated pass; of the resulting
41 candidate claims, every one was manually inspected against the arm-specific input it
actually describes. This is a claim-level audit, not line-by-line prose reading of all
58 files -- prose that makes no `param = value`-shaped assertion was not checked.

## Findings: parameter-claim table

| run dir | claim | file:line | actual | verdict |
| --- | --- | --- | --- | --- |
| `bench/cs_mhd_prod3_wb` | NOTES.md: "removes the WB difference only: `wellbalance_dynamic=true, wb_x1=true, wb_option=polytropic, wb_cache_every=10` (sp_mhd_prod3's block)" -- implies the run reproduces sp_mhd_prod3's WB block in full | `cs_mhd_prod3_wb/NOTES.md:5` | `deep_hot_jupiter.athinput` sets the 4 named keys (lines 181-184) but **does not set `cs_wellbalanced_src`** anywhere outside a comment; it defaults to `false` (`src/coordinates/coordinates.cpp:43`) | **CONTRADICTION** (this is the task's cited known example, now precisely located: the run's own header comment, still copy-pasted from `cs_mhd_prod3`, says WB is "OFF" and lists `cs_wellbalanced_src` among the omitted keys at lines 6-7/156-158, while the body sets the other 4 WB keys ON -- an internally inconsistent input whose comment and settings disagree, and whose NOTES.md then describes it as reproducing sp's *full* WB block when the cubed-sphere geometric WB source term never turns on) |
| `bench/cs_wb_cache/every10` and `every1` | (self-correcting) README.md:43-58 explicitly states that `cs_mhd_prod3_wb/NOTES.md` (as read in an earlier session) claimed a fifth setting, `cs_wellbalanced_src = true`, and verifies the input has no such key | `cs_wb_cache/README.md:50-56` | `cs_wb_cache/every10/deep_hot_jupiter.athinput`: `cs_wellbalanced_src` appears only in a comment (line 51 region), confirmed absent as a live key | **CONFIRMED, already documented** -- `cs_wb_cache/README.md` itself is the correction; it matches the current `cs_mhd_prod3_wb/deep_hot_jupiter.athinput` state exactly. (Note: the *current* text of `cs_mhd_prod3_wb/NOTES.md` on disk today only lists 4 keys, not 5 -- the 5th-key claim `cs_wb_cache` is correcting may be from an intermediate edit of that NOTES.md, or from verbal/session context rather than the file as it reads now. Either way the underlying fact -- `cs_wellbalanced_src` was never live in that arm -- is real and worth carrying forward into any new WB-on cs production input, exactly as flagged in the `deep_hot_jupiter_cs_prod4.athinput` header written for Task A.) |
| `bench/cs_mhd_prod3/NOTES.md` | "`problem/rt_ck` ... is already `true` in cs_mhd_prod2's input and is kept. It forces `rt_split = true`" | `cs_mhd_prod3/NOTES.md:75-76` | `deep_hot_jupiter.athinput` line has `rt_split = false`; `src/pgen/deep_hot_jupiter_rt.cpp:529-532` does force `rt_split = true` at startup whenever `rt_ck && !rt_split` | **NOT a contradiction** -- correctly describes a runtime override of an input-file default; the raw input value and the described *effective* value are different things by design. |
| `bench/rg_box/NOTES.md` | "in the top `cool_depth = 0.3 H_p`" | `rg_box/NOTES.md:170` | `rgbox.athinput:117`: `cool_depth = 1.866568e10   # 0.3 H_p` | **MATCH** -- 1.866568e10 is exactly `0.3 * H_p` per the input's own comment. |
| `bench/sp_excess/mhd_thuni`, `hyd_thuni` | "`f_stretch_theta = 1.0e-3`" | `sp_excess/NOTES.md:85-91` | `mhd_thuni/deep_hot_jupiter.athinput:122` and `hyd_thuni/deep_hot_jupiter.athinput:122` both read `f_stretch_theta = 1.0e-3` | **MATCH** (only a sibling arm, `mhd_lowbeta`, keeps the production `3.0` -- comparing the claim against that sibling is what makes the naive automated pass misfire). |
| `bench/cs_deep_ablate/{nowb,nogs07,etacap,etalow}` | per-arm one-line changes: `cs_wellbalanced_src = false`, `cs_gs07_emf = false`, `ohmic_resistivity = constant` + `eta_ohm_const = 1e13`, `max_eta = 1e12` | `cs_deep_ablate/README.md:4,5,21,24` | all four verified present verbatim in the correct arm directories (`nowb/`, `nogs07/`, `etacap/`, `etalow/deep_hot_jupiter.athinput`) | **MATCH**, all four. |
| `bench/ck_sph_ab` | `ck_spherical=false` / `true`/`false` / `true`/`true` across `off`, `sph`, `sphbeam` | `ck_sph_ab/README.md:97,172,179` | `off/`, `sph/`, `sphbeam/deep_hot_jupiter.athinput` each carry the matching `ck_spherical`/`ck_beam_sph` pair | **MATCH** per-arm (the file's own README already documents, and this audit confirms independently, that `off`'s *result* is not a valid control for unrelated reasons -- see `tests_gpu_ab_0922/README.md` -- but the **input values themselves** are exactly as claimed). |
| `bench/cs_deep_ablate` | "same binary as `cs_prod_mhd_rot` (md5 97bf06ab...)" | `cs_deep_ablate/README.md:3` | `md5sum cs_deep_ablate/athena` = `97bf06abf0147922b986dd237fa230e4` | **MATCH** |
| `bench/cs_ablate_r` | "ORIGINAL GPU binary of cs_hyd_rs (commit f6f0d6f4, md5 b40447e6dad9409d674f7fccbdd5290e)" | `cs_ablate_r/NOTES.md:22` | `md5sum cs_ablate_r/{old_cache1,h_old_cache1}/athena` = `b40447e6dad9409d674f7fccbdd5290e` | **MATCH** |
| `bench/sp_mhd_prod`, `sp_mhd_prod3`, `cs_n64_fixed` | `md5 ca06b01c`, `md5 4fa129a1`, `md5 efd92190` in `BUILD_COMMIT.txt` | see each `BUILD_COMMIT.txt` | `md5sum` of each dir's `athena` matches to the stated 8-hex-digit prefix | **MATCH**, all three |
| `bench/bstar_fecz` | "`bc_mode = 2` at both x1 walls, not smoke/'s `bc_mode = 0`" | `bstar_fecz/README.md:358` | no `.athinput` under the directories walked (depth<=3) sets `bc_mode = 2`; closest matches (`final_gate/A1{b,c}.athinput`) set `bc_mode = 3`, and `ic/probe_rt.athinput` sets `bc_mode = 0` | **INCONCLUSIVE** -- the specific run this sentence describes was not identified among the `.athinput` files present (it may be an unlabeled/overwritten arm, or a run whose input has since been deleted/replaced); flagging rather than calling it a contradiction, since the referenced file may simply be gone. |

## BUILD_COMMIT.txt commit-hash sanity check

Every `BUILD_COMMIT.txt` leads with a **git commit hash** (7-8 hex chars), not a binary
md5; comparing that value against `md5sum <dir>/athena` (as a first, naive pass did) is a
type error and always "mismatches" -- that is not a finding. Only the *explicit* `md5 …`
substrings inside these files are genuine binary-identity claims, and every one of those
found in the 58-file corpus (5 total, listed above) checks out against the binary that is
still on disk. No BUILD_COMMIT.txt / binary contradiction was found. Several
`BUILD_COMMIT.txt` files reference binaries that are no longer present in the directory
(e.g. `cs_mhd_prod3_fofc`, `ck_mhd_b3` retains only a git-commit BUILD_COMMIT.txt with no
persisted md5 to check) -- these are unverifiable, not contradicted.

## Summary

- 58 NOTES.md/README.md files audited (100% of the corpus at depth<=2 under `bench/`).
- 41 automatically-flagged parameter claims, all manually resolved.
- 1 confirmed real contradiction: **`bench/cs_mhd_prod3_wb`** -- its own header comment
  and NOTES.md describe reproducing sp_mhd_prod3's well-balanced block, but
  `cs_wellbalanced_src` (the cubed-sphere-specific piece of that block) was never set in
  the actual input and defaults to `false`. This is already independently caught and
  documented by `bench/cs_wb_cache/README.md`.
- 0 contradictions in the 5 binary-md5 claims checked (all matched `md5sum` exactly).
- 1 inconclusive claim (`bstar_fecz` `bc_mode = 2`) whose specific run input could not be
  located.
- All other flagged claims were confirmed correct once matched to the arm directory they
  actually describe.
