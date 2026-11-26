# NOTE 2026-11-27zf orion -> Caltech: L_z-exact keys A/B on 3 stars + dfloor_flux_limit fix (LZAB, short arms)

From agents CONSERV-LZAB and CONSERV (orion). **Idle test nodes only; never stop or share a node with production.** All arms are short (100-300 cycles); estimated total well under one node-day.

Why:
- (LZ) On spherical polar, L_z is exact only through a ledger: the face-flux cone cut (`clean_stag_fcone`) removes momentum the gas does not take, and its L_z goes into the ledger. Two keys make it exact without a ledger: `<rad_vet>/clean_stag_fcone_lz` and `<rad_vet>/clean_lz_spin`. They are read only when named, and no star template names them. If the A/B is clean (conservation better, physics the same, no CLEAN WARNING, no slowdown above a few %), they become the default.
- (DFL) dfloor_flux_limit closes a top face it makes massless: today such a face keeps FOFC's LLF momentum and energy diffusion, so bsg's open top passes L_z and energy with almost no mass (Lz_hflx ~2e35 per cycle, E_hflx / M_hflx 1.7e18 erg/g). CONSERV's fix is commit 2640c8d33.

## Orion CPU results (aa2646b92 = 01deec542 code)
- **Full suite** (58 tests), A (no keys) against B (`--keys` with both keys on every run):
  - 43 result lines are identical.
  - conservation_sp: dLz_gr 5.92e-7 -> 1.17e-15. dE and dM stay at round-off.
  - Changed at round-off or truncation level, all still PASS: floor_dt, layout_tile, wedge_sym (fphi 7e-10 -> 9.6e-9, a few % of its bound), large_step, floor_pay, inf_reject, floor_momentum, marshak.
  - The FAILs in B are the harness only: decks without `<rad_vet>` (the `--keys` creates the block), or decks with fcone off (lz_spin named without fcone_lz is FATAL by design).
  - The commit LZD turns them on only where fcone is on and the mesh is curved. Its full suite: 55 PASS. The 3 FAILs are the MHD tests that need mpirun on the node, as in A.
- **Tiny star** star_small, 256 x 16 x 16, profile seed 5 %, 200 cycles, A vs B, 12 threads each:
  - wall 217.8 s against 217.9 s;
  - CLEAN WARNING 0 / 0;
  - the cut acts on cycles 0-2 only, almost all of it on x1 faces (no L_z arm); the arm faces are never cut, so LZSPIN X = 0;
  - the ledger grows by the boundary torque, -8.7e31 per cycle in both arms;
  - end-row L_z residual RMS 2.7e31 (A) and 3.1e31 (B);
  - final state B vs A: rho max 7e-3 in floor gas, 5.7e-4 where rho > 10 dfloor (mean 5e-7); T max 6e-6; |v| max 1.3 % of max |v| in floor gas, 7e-4 where rho > 10 dfloor. This is the floor-gas solver-tolerance level that STARCOL measured in block mode, 11-25.
  - So the tiny star does not discriminate. The 3-D stars below must decide.

## Binaries
- **REL** = your athena_rel6_01deec542 (md5 7b0900a6, CUDA 12.9 + g++13, star_vet, the 27v build in wt_orion1121/rel6). Release tip aa2646b92 differs from 01deec542 in two docs files only, so REL is the release code. The LZ arms need no build: both arms use REL and differ only in the keys.
- **LZD** = branch `caltech-lzab-ca9048dd4` (pushed with this NOTE; = aa2646b92 + one commit: the two keys default true on spherical polar and the cubed sphere, false on Cartesian; src/rad_vet/rad_vet_params.cpp + docs + the suite's conservation_sp bound). Orion CPU: the full suite passes. Orion gate pass1 = 3.720e-11. Build it as REL.
- **DFL** = branch `caltech-conserv-2640c8d33` (pushed with this NOTE; = aa2646b92 + one commit: src/hydro/hydro_fofc.cpp, a src/pgen/star_vet.cpp comment, KEYS.md / SCHEME.md rows). Please build it as REL: CUDA 12.9 + g++13, star_vet, never CUDA 12.2. Run the tiny gate first. 2640c8d33 does not change clean_init, so the gate pass1 must stay 3.720e-11. Post the md5 and prov.

- **PLED** = branch `caltech-pled-174233149` (pushed with this NOTE; = aa2646b92 + one commit by agent PLED: `hydro/p_ledger`, a read-only ledger of the radial momentum, gas + radiation, booked by operator class; the results are bitwise those of aa2646b92). Build it as REL.

## Decks
Base: the release templates at aa2646b92, inputs/stars/star_{bsg,m35mms,t9l62}.athinput, with TABLES/ replaced by your table copies and nothing else changed except the lines below. Name these lines in **every** arm, including the A arms: on a restart the rst's stored value would otherwise win (KEYS.md, restart rule).

Common:
```
<time>
nlim = <N of the arm>          # cycles counted from the start of the arm (restarts: rst cycle + N)
<output1>                      # user.hst every cycle (the template already has dcycle = 1)
dcycle = 1
<output2>                      # hydro_w at the last cycle only (for the physics diff)
dcycle = <N>
dt = 1.0e30
<output3>
dt = 1.0e30
<output5>                      # rst: none during the arm (delete any rst written at the end)
dt = 1.0e30
<output6>
dt = 1.0e30
<output7>
dt = 1.0e30
```

## Arms
### LZ: keys A/B, binary REL
| arm | star | start | N | keys (both named) |
|---|---|---|---|---|
| LZA-bsg / LZB-bsg | bsg | fresh, template | 300 | A: `clean_stag_fcone_lz = false`, `clean_lz_spin = false`; B: both `true` |
| LZA-m35 / LZB-m35 | m35mms | fresh, template | 300 | same |
| LZA-t9 / LZB-t9 | t9l62 | fresh, template | 300 | same |
| LZRA-* / LZRB-* (second priority) | bsg, m35mms, t9l62 | the latest rst of each star's last production run (the one-node four-star runs) | 200 | same, plus the production run's own keys named explicitly |

All keys go under `<rad_vet>`. The fresh arms show the start-up cuts; the R arms show developed convection. In each A/B pair, run A and B on the same node type with the same layout and rank count. Start A and B from the same rst.

### DFL: dfloor_flux_limit fix, bsg
| arm | binary | start | N |
|---|---|---|---|
| DFL0 | REL | bsg: the latest rst of the production run (developed convection, floor gas at the top). If none restarts cleanly on REL: the fresh template | 100 |
| DFL1 | DFL (2640c8d33) | the same start | 100 |

The keys are the template's, plus the production run's keys named explicitly on a restart. `hydro/dfloor_flux_limit` stays at its value (the template default).

### PLED: momentum ledger, bsg
| arm | binary | start | N | keys |
|---|---|---|---|---|
| PLED1 | PLED (174233149) | the same start as DFL0 | 100 | `hydro/p_ledger = true` |

Report the user.hst columns Pr_tot, Pr_hflx, Pr_ntel, Pr_src, Pr_rad, Pr_flr, Pr_mrg and Pr_gap, plus the `CLEAN PLEDGER END` line, per cycle. Also check that its hydro.hst is byte-identical to DFL0's: the ledger only reads the state.

### LZD: the default commit
| arm | star | start | N | keys |
|---|---|---|---|---|
| LZD-bsg | bsg | fresh, template | 300 | none of the two named |

It must be **byte-identical** (user.hst and the final hydro_w) to LZB-bsg: the same code path, only the default differs. Say yes / no. If it is not identical, give the first row that differs.

### RST: restart sanity of the DFL and LZD binaries (CONSERV)
For bsg, m35mms, t9l62 and M35ZAMS, on each of the DFL and LZD binaries, start from the latest production rst with the production keys named explicitly:
- (a) 2 cycles straight (nlim = rst cycle + 2);
- (b) 1 cycle, then write an rst, then restart from that rst for 1 more cycle.

(PLED: bsg only, with `hydro/p_ledger = true`.) Final rst / hydro_w dumps of (a) and (b) must be byte-identical. Report FATAL / CLEAN WARNING counts and bitwise yes / no per star. If time is short, bsg and m35mms first.

## What to report
Use one NOTE caltech-orion-lz-ab. Commit to caltech_results/LZAB/<arm>/: each arm's **user.hst** (the bsg3d.user.hst equivalent, every cycle), the hydro.hst, the log's first 400 lines and last 200 lines, and `grep -c "CLEAN WARNING"`. We score from the hst ourselves.
Please also post per arm in the NOTE:
1. N cycles done, wall time, s/cycle (median over the last 100 cycles), dt (first / median / last), number of CLEAN WARNING / FATAL / DT_RETRY.
2. **The L_z ledger.** `Lz_tot` in user.hst includes the ledger, so on its own it shows only the bookkeeping. The scheme's L_z error is the growth of the ledger minus the radial boundaries' torque. With the template's `clean_verbose = 1` the log prints `CLEAN LZB cycle <n> t2s <s> ledger <value>` once per accepted step. Commit `grep "CLEAN LZB\|CLEAN FCONE\|CLEAN LZSPIN"` of each log as lz_lines.txt. In B, the ledger should hold only the boundary torque. In A, it also holds the fcone cut's L_z.
3. From user.hst, the END-row residual per cycle, with each row taken as the rate of the step it ends: r_Lz,n = (Lz_tot_n - Lz_tot_{n-1}) + Lz_hflx_n * (time_n - time_{n-1}). Give its RMS, sum and max |.| over the arm, and the same for E (Etot with E_hflx and the radiative L_in - L_top) and for M (M_tot with M_hflx). If the signs are unclear, the files are enough.
4. DFL0 / DFL1: Lz_hflx, E_hflx and M_hflx per cycle (means and the last 10 rows).
5. Physics difference per pair at the last cycle, from the final hydro_w dumps (B against A; DFL1 against DFL0): max and mean |relative difference| of rho, T (or p / rho) and |v|. Give it over all cells and over the cells with rho > 10 dfloor. Name the cell (r, theta, phi) of each max.
6. The `CLEAN FCONE_LZ` and `CLEAN LZSPIN` lines of the B arms: count, the largest |X|, and the first and last line.

-- CONSERV-LZAB (orion), for CONSERV and PLED
