# NOTE vimpbcdef (worker R, 2026-10-12): `implicit_vimp_bc = true` as the fresh-run default

Work dir /raven/ptmp/jinma/vimpbcdef_1011.
Branch `vimpbcdef-1011` (worktree /raven/u/jinma/ATHENAK/wt-vimpbcdef), from prodfinal-1011 422c2af7. Not pushed.

## 1. Change (commit 67b46926)

- `src/rad_m1/rad_m1_implicit.cpp` (the key's parse, just after `implicit_bc_advect`):
  - **implicit_vimp on:** `impl_vimp_bc = GetOrAddBoolean("rad_m1","implicit_vimp_bc", !restart_run)`.
    - A fresh run takes true.
    - A restart whose file lacks the key keeps the old value false.
    - The resolved value is echoed into the parameter set, so later restarts keep it.
    - Explicit input wins.
  - **implicit_vimp off:** unchanged. The key is read only when named, and it is not echoed (it acts nowhere without
    vimp).
- The key still acts only with `implicit_bc_advect`, at a Marshak x1 end face (`vbc = impl_vimp_bc && impl_bc_advect`,
  ImplicitVimpBuild).
- Comments are updated in `rad_m1_implicit.cpp` (the kernel note) and `rad_m1.hpp`. There is no other code change.
- `src/pgen/he_star_m1.cpp` is untouched.

Binaries (bin/), built by build.sh (job 31088314):
- athena_cpu_67b46926, md5 27e2c9c6a52e64f3a9c6017be502fdb0;
- athena_he_cpu_67b46926, md5 4331a343be381c3bc30d48ae1c150395.

The "old" binary is /raven/ptmp/jinma/consfix_1011/bin/athena_{cpu,he_cpu}_1fbe8621. Its code equals 422c2af7
(422c2af7 is docs only).

## 2. Checks

Note on method: AthenaK .bin files carry the input-parameter text in their header. Two decks that differ only in a
named key therefore never `cmp` equal. "Bitwise" below means every data array of every dump (bin) or every byte of every
.tab is identical (ana/bcmp.py).

### 2.1 Key named vs default (rlat_T, rswv_T, rswv_cart; level 1)

| pair | result |
|---|---|
| old + `implicit_vimp_bc=true` vs new + `implicit_vimp_bc=true` | **bitwise** (3/3 cases) |
| new, no key vs new + `implicit_vimp_bc=true` | **bitwise** (3/3) |
| new + `implicit_vimp_bc=false` vs old, no key | **bitwise** (3/3) |

### 2.2 Restart (rlat_T, rst at t = 2, run to t = 4; rst.sbatch, job 31088438)

| run | result |
|---|---|
| new binary restarted from the OLD rst (no key) vs old binary restarted from it | **bitwise** (false kept) |
| old-rst restarts vs the old fresh run | bitwise (the restart is exact) |
| new fresh rst header | `implicit_vimp_bc = 1  # Default value added at run time` |
| new from the NEW rst vs old from the NEW rst vs the new fresh run | **bitwise** (true kept, read by both) |
| new from the old rst vs new from the new rst | differs (8e-2 max rel in v), as it should |

### 2.3 Where the default flips (fresh, no key; old = false, new = true)

T self-convergence ladders, last dump (h.py tord, sptestsC norms). These are vimp_bc alone, with no `time2_start=slope`.

| case | old default (orders 1-2 2-4 [4-8]) | new default |
|---|---|---|
| rlat_T (lateral wave, Marshak top) | dens 1.54 1.34; velx 1.54 1.33; E 1.54 1.33; F1 1.54 1.33. L1(E, 4-8) 7.5e-5 | **dens 1.95 1.98; velx 1.92 2.04; velz 1.93 2.04; E 1.95 2.10; F1 1.93 2.07; F3 1.93 2.03. L1(E, 4-8) 2.2e-6** |
| rswv_T (v_r 0.01) | dens 1.11 1.02; E 1.08 0.97; F1 1.07 0.95. L1(E, 4-8) 6.6e-5 | dens 1.21 1.05; E 1.48 0.66; F1 1.40 0.65. L1(E, 4-8) 9.2e-6 (7x lower; the known rswv floor, H section 3.3, open) |
| rsw_T (static pulse) | dens 2.02 1.90; E 1.76 1.41 | dens 2.03 1.91; E 1.77 1.41 (equal; max diff 7e-9 rel) |
| rswv_cart base (no flow) | 1.77 1.41 1.01 (E, levels 1-16) | identical digits (diff 1e-11) |

- The old-default numbers reproduce H's H0 (prodfinal default) rows exactly. The data match H0 to ~1e-8 rel in v: the
  code moved from 18ff91c8 to 422c2af7.
- rswv_T new = H's "slope + vimp_bc" E 9.2e-6.
- rlat_T with `+time2_start=slope` on the new binary (NS, levels 1-8) gives dens 1.94 2.16, E 1.94 2.26, F3 1.89 2.07.
  These equal H's FIN table (the data differ by ~1e-9 rel from FIN, the same code drift).
- **New finding.** On rlat, vimp_bc alone already brings 1.5/1.3 to 1.95/2.0. slope adds little there: L1 1.7e-6 vs
  2.2e-6 at 4-8.
- Irradiated cloud (suiteC `cloud REC r12g2`, Marshak ends, vimp on, so the default flips): not bitwise, but at round-off
  (F1 1.2e-10 abs of 7.3; v 5e-13). This matches J's AG Car outflow finding.

### 2.4 Must not change (old vs new, no key)

| run | why | result |
|---|---|---|
| rw_p10t100, rw_p10t1 (P = 10 waves, vimp on) | periodic, no Marshak face | **bitwise** |
| eq3t REC (moving equilibrium, vimp on) | periodic | **bitwise** (402 dumps) |
| radiation shock M 1.2 REC n512 (vimp on, Marshak both ends) | deck sets `implicit_bc_advect=false` | **bitwise** |
| Marshak REC c dt/dx 1 and 10 | `implicit_vimp=false` | **bitwise** |
| AG Car A reduced, FRESH deck A0, 60 cycles, he binary (agc/agc.sh, job 31088551) | top cells inflow, E_bath = 0 (H, J) | **bitwise**: old default = old + key true = new default; hst md5 equal; Picard 3.02/6, 0 NONCONV |

Cost: unchanged. Run times are equal within noise (cloud 566 vs 598 s, rs 319 vs 306 s).

## 3. Production

- Fresh AG Car decks get true by default. prod/agcar_B_fresh_v2 already names it; A/B fresh decks without it now
  resolve true (bitwise on AG Car, section 2.4).
- Restart chains (restart_keys_B_runB.txt does not name it) keep false, as the rule intends. Naming
  `implicit_vimp_bc=true` on a restart command line would turn it on: a decision for the coordinator.

## 4. Files

- ana/:
  - r.sh (h.py runner, O/N binary);
  - s.sh + run1.sh (suiteC runner copies);
  - h.py (H's harness, W changed);
  - bcmp.py (data compare, header excluded);
  - cmpdir.sh (raw byte compare; it misleads on bin headers).
- Job lists: jobs_A{O,N}.txt and jobs_B{O,N}.txt (pool.sbatch / poolS.sbatch, jobs 31088315/410/316/437).
- Runs: runs/<case>/{Od,Nd,Ot,Nt,Nf,NS}_<lev>, runs/rlat_Tr2/{Od_1,Nd_1,ROo,RNo,RNn,ROn},
  runs/{cloud,rs_m1.2,marsh,eq3t}/REC_{O,N}*, agc/runs/j.31088551/{Od,Ot,Nd}.
- base/: copies of H's decks; rlat_Tr2 = rlat_T + rst every t = 2.
