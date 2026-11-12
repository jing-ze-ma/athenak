# TASK UPDATE 2026-11-12 viper -> Orion: idortd6 (1d5716074) runs on the d6fast-1112 build, NOT 421a3d024

The user changed the gated object: it is now **d6fast** = the D6 keys on the build of `fork/d6fast-1112` at commit
**5f06e171b236a2215b24fbcdfacab21b927039a5** (= fork/fast2nd-1112 e1ce0a1e1 on 421a3d024 + the fast vet_gd shell kernel 79030febf
+ `#pragma clang fp contract(on)` 58b9ca152; the source tree is identical to fork/fpon-d376c311c 01516be4a).
Everything in TASK-2026-11-12-orion-idortd6.md (1d5716074) stays the same (lists, keys, decks, drivers, analyzers, pass rules,
return NOTE) except the binary:

1. Build: a NEW worktree, never touch main and do no destructive git in existing worktrees:
   ```
   cd /orion/ptmp/jinma/Athenak && git fetch fork && git worktree add --detach wt_d6f 5f06e171b236a2215b24fbcdfacab21b927039a5
   ```
   Configure `wt_d6f/build_none` with exactly the CMake options of `wt_fin/build_none` (copy them from its CMakeCache.txt: same
   compiler, MPI, OpenMP, build type, no PROBLEM), build, and record the md5 in the NOTE.
2. Run both parts of 1d5716074 with `RW_EXE=/orion/ptmp/jinma/Athenak/wt_d6f/build_none/src/athena` (part 1) and
   `EXE=` the same binary (part 2). Use NEW run dirs with a `d6f` tag (`runs_finD6f_<c>`, `fin_D6f_<c>.out`, `$O/curv/<cell>_D6f`,
   `ev_all_d6f.txt`) so nothing of a 421a3d024 run is overwritten. Smoke first as in the TASK.
3. If you have already started or finished 1d5716074 on 421a3d024: keep those results (report them as "D6 on 421a3d024"),
   and rerun everything on the d6fast build. The verdict that counts is d6fast's.
4. The FIN references stay your FIN runs on the 421a3d024 build (wt_fin). Also give d6fast vs your 421a3d024 D6 if you have it
   (expected: identical or round-off; on gcc the new shell kernel is a different code path).

Return: in the NOTE of 1d5716074 (NOTE-2026-11-1x-orion-viper-idortd6.md), with a clear d6fast column and PASS / FAIL per part.
