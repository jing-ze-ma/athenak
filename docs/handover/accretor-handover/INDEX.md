# Accretor (ry_per_accretor / Plaskett) hand-over from viper, 2026-10-07

viper ran this project 10-06 .. 10-07; from 10-07 it moves to DeltaAI. Everything here is a verbatim copy
of viper's working directory `/viper/ptmp2/jinma/accretor_1006` (small text files only), plus this index.
Numbers quoted below come from the copied notes or the named run logs; where I add a derived number it says
so.

## Files in this directory

| file / dir | what it is |
|---|---|
| DESIGN.md | 10-06 design for RY Per (2-D r-phi, rotating frame about the accretor, absorbing stage-1 surface), literature parameters, ballistic stream numbers, grid, measured quantities, staged plan, open decisions |
| IMPL.md | stage 1 (absorbing rigid surface) implementation + tests: at-rest atmosphere, cold ballistic stream vs ballistic orbit, spin 1.0 vs 7.2, GPU smoke + throughput |
| IMPL2.md | stage 2 (resolved polytropic envelope, `problem/inner = envelope`): trade-off, what was built, failed attempts, ambient/floor redesign, all RY Per GPU tests and collapse diagnoses, hot-star test, spin pair, option 3 (resolved real photosphere) |
| CANDIDATES.md | literature scoping of disc-forming Algols and massive gainers (phi Per progenitor, Plaskett progenitor, observed systems), cost, ranking |
| PLASKETT.md | Plaskett's-star progenitor setup (Wade+2026), env_rho_ph = 0.3 estimate, options 1/2, radiation plan, q10/q11 runs and collapse diagnosis |
| PENDING_CMDS.txt | the exact commands for qd1/qd2, q11, q12 (written when blocked by the permission check; all later run) |
| msg1-3.txt | short messages written to the coordinator on 10-06 |
| build.sh, build_dbg.sh, patch_iso.py | build scripts (GPU rocm 7.2 / CPU debug) and a one-off patch helper |
| bin/*.athinput | every input used: ry_per_accretor_{f5dfc1ef,1e8d0ff8,env1..env8} (RY Per), plaskett_env9..12, 12f |
| keys/<tree>_<arm>.txt | the command-line key file of every arm (gpu2/*, gpu/*, tests/*, tests2/*) |
| jobs/*.sbatch | the submit scripts (B= binary, INP= input; arms are the run directories given as arguments) |
| scripts/ | analysis + design: roche_stream.py (ballistic L1 stream), grid_design.py (+ grid_plaskett_500.txt), ana_env.py, budget_env.py, dtcell.py, dtcell2.py, mkdiag.py, mkfofc.py, mkhot*.py, mkplaskett*.py, mkres.py, prof.py, tser.py, col.py, where.py, env_tradeoff.py |
| tests/ | stage-1 CPU analysis (ana.py, budget.py) and run scripts (bitwise_stage1.sh, run_cpu.sh, run_env_cpu.sh) |
| page/ | mkfigs.py, mkfigs2.py (figures; the Plaskett page script mkfigs_plaskett.py is already in plaskett-1007/page on this branch) |
| cand/ | candidate scoping scripts and their outputs (cost.py, cost2.py, scan.py, ls75.py, epochs, systems tables, stream_runs.txt) |
| lit_files.txt | file list of viper's `accretor_1006/lit/` (arXiv PDFs + pypdf text dumps the notes cite by line number) |

Not copied (big, viper only; DeltaAI cannot pull them): binaries
`/viper/ptmp2/jinma/accretor_1006/bin/athena_ryper_{cpu,gpu72}_*`, build trees `build_cpu`, `build_cpu_dbg`,
`build_gpu72`, `lit/*.pdf`/`*.txt`, all run output under `gpu/`, `gpu2/`, `tests/`, `tests2/` (bin dumps,
rst, hst; e.g. gpu2/q12_s45 1.4 GB, q12_a1 2.4 GB, q10_s1 0.95 GB, t7_s1 1.2 GB), figures `page/*.png`,
`page/plaskett_test/`.

## Science goal and quantities (DESIGN.md "Measured quantities")

How much angular momentum the stream gives a mass-gaining star, and whether a fast-spinning gainer still
spins up (RY Per: observed Omega_* = 7.2 Omega_orb; the user later moved the target to massive gainers,
the Plaskett progenitor). Measured through `user_hist` (22 columns, IMPL2.md):
- Mdot_in (stream window), Mdot_acc (through r_in or R_acc), Mdot_out;
- Jdot through the star, inertial, about the accretor, plus the surface stress (Jstr);
- **j_acc = Jdot/Mdot_acc vs j_Kep(R_acc)**, spin-up time J_*/Jdot;
- disc vs direct impact: phi distribution of Mdot/Jdot at r_in, Sigma(r), v_phi/v_K.
Results are per unit Mdot as long as the run is scale-free (stage 1); the envelope mode is not scale-free
(env_rho_ph, see NOTE section 2).

## Code history (branch accretor-1006)

| commit | change | why |
|---|---|---|
| f5dfc1ef | pgen ry_per_accretor: L1 stream on an absorbing rotating star (stage 1) | DESIGN.md stage 1 |
| (1e8d0ff8 binary) | stage-1 binary after the driver NaN-scan fix | IMPL.md |
| fe16e65a | stage-2 resolved-envelope mode `problem/inner = envelope` (n = 3 polytrope in Phi_wb, ideal gas, etotgrav + wellbalance_dynamic + wb_x1, thermal relaxation, hot hydrostatic ambient, floor sponge); binaries env1..env5 | stage 1 cannot measure a spin-dependent torque (IMPL.md (c)) |
| 6588011b | keys env_cs_stream, env_r_hot (defaults reproduce fe16e65a); binary env6 (md5 77a17639...) | hot-star test: stream injected cold (15.5/21 km/s) while the star is hot |
| ea04cd19 | key env_r_spin (the star's rotation and the spin term of Phi_wb end there); binary env7 (md5 a5ad0a18c1d3a73e4cb41e833f1821d5) | t7_s72: at spin 7.2 the spin-inclusive potential peaks near corotation, so the hot atmosphere co-rotated to r_out and was flung out |

All Plaskett runs (q10, q11, q12) use binary env7.

## Run directories on viper (all under /viper/ptmp2/jinma/accretor_1006)

### RY Per, stage 1 (absorbing surface; IMPL.md)
| run | what | outcome / lesson |
|---|---|---|
| tests/a1w*_hot150, a1w2_hot300, a1_hot150/300 (CPU) | at-rest atmosphere, wall vs diode | closed wall at rest to 2.7e-4 c_s (300 km/s); absorbing diode drains a 150 km/s atmosphere (expected) |
| gpu/a2_rin406, gpu/a2_rin2 (jobs 12105359/8) | cold (1 km/s) ballistic stream | ridge within 0.45 deg of ballistic; impact phi 74.09 vs 73.65 deg; r_min 2.774 vs 2.695 (+3 %) |
| tests/c_spin1, c_spin72 (CPU, quarter res) | spin 1.0 vs 7.2, stream on | j_acc/j_Kep = 1.0647 in both: **an absorbing surface makes the torque spin-independent** -> stage 2 |
| gpu/smoke, smoke2, gpu2/smoke, gpu2/tput | smoke + throughput | 3.41-3.44e8 zone-cycles/s/node; 23 min/orbit/node |
| tests/dbg*, tests/full_cpu | CPU debug (bounds check) | no OOB in 210 cycles |

### RY Per, stage 2 (envelope; IMPL2.md sections 3-6; binaries env1..env7)
| run | what | outcome / lesson |
|---|---|---|
| tests2/dbg_a..dbg_h3, a_env_*, a2_env_*, ah_env_*, amb_h1, smoke_env, bw_old/new (CPU) | failed attempts | WB on the unsupported ambient (H_p << dr) -> collapse cycle 5; WB only below R_acc with analytic IC -> NaN; top-anchored discrete HSE -> unbalanced; quarter res too coarse; ah_env dt collapse -> ambient redesign (hot hydrostatic ambient, env_amb_rho, env_amb_k, sponge) |
| gpu2/a2_*, a3_* | log-rho interpolated relaxation target; ambient 1e-4 | evaporation wind 1000-2900 km/s; ambient P interacts with the cold top |
| gpu2/g_half, g_halfslow | two-phase T attractor | precursor jet 1700-1800 km/s |
| gpu2/f_half, f_full, h5_inst | instant relaxation (1e-9) | 375x pressure jump, collapse (f_half cycle 1032, h5_inst cycle 99) |
| gpu2/h5b_str | env_r_top = R_acc | collapse t = 0.006 |
| gpu2/e_amb6, e_both, e_tmid, d3_s1, dg_s1, a9/a10_s1.0 | envelope-only variants at half res | (no separate note; keys in keys/) dg_s1 FATAL on an unknown key time/dt_min (binary predated it) |
| gpu2/c_s1.0/7.2, csm_*, c4_*, b4*, s9/s10 | stream runs on early ambients | c_s1.0 collapse t = 0.057; c4 collapse t = 0.0119; s10 collapse t = 0.0697 |
| gpu2/h5_a1, a4_s1.0/7.2 | envelope alone 0.2 orbit (test a) | PASS; s 7.2 Jenv drift -2.35e-4 (tidal torque), wall-cell feature 5 km/s unexplained |
| gpu2/h5_str, h5_imp | stream, WB everywhere | collapse t = 0.0499 at the photosphere |
| gpu2/i_base, iA-iD (jobs 12108343/95) | restart diagnosis | deterministic; **WB reconstruction in cold under-resolved atmosphere cells (H_p = 0.13 dr) under non-hydrostatic inflow** dumps energy into the plain cell above; wb_rmax = R_acc survives; faster relaxation collapses earlier |
| gpu2/a6_s1.0/7.2, s6_s1.0, s7_base, s7_r | wb_rmax = R_acc | envelope OK; stream passes first impact, collapses t = 0.069: one-step positivity loss at the WB/plain face under Mach ~10 shear |
| gpu2/t6sm, t6_a1, t6_s1 (env6) | hot star c_ph 85.3 km/s (H_p = 4 cells), WB everywhere | first run through the impact (t = 0.31); impact phi 76.5-81.4 deg; ~90 % of stream AM piles up outside R_acc; not steady |
| gpu2/t7_s1, t7_s72 (jobs 12108898/9) | hot-star spin pair to t = 0.511 | s1: j of accreted mass 0.67 j_K (window-to-window 0.07-1.19); s7.2 INVALID (corotation fling) -> ea04cd19 env_r_spin; rerun on hold |
| gpu2/r8sm, r8_s1, r8_a1 (env8) | option 3, resolved real 15.5 km/s photosphere, nx1 900 | smoke clean (job 12108943, dt 4.5e-7, 2.50e8 zc/s/node); production cancelled (12109190/1, too expensive: ~15.5 h/orbit/node) |

### Plaskett progenitor (PLASKETT.md; binary env7)
| run | job | what | outcome / lesson |
|---|---|---|---|
| (plaskett_env9) | — | option 1, hot surface c_ph 64.8 | never run |
| gpu2/q10sm, q10sm_ns | 12109320 | option 2 smoke (500 x 4 x 2048, real surface c_ph 20 km/s, WB everywhere) | clean; 2.43e8 zc/s/node |
| gpu2/q10_a1 | 12109350 | envelope only 0.2 orbit | PASS: envelope max dv 0.024 km/s, Menv -8.3e-6, Jenv -1.6e-5 |
| gpu2/q10_s1 | 12109351 | stream, 2 apu nodes | **dt collapse t = 0.03531 (cycle 19177, 0.077 orbit)** |
| gpu2/qd1, qd2 | 12109481 | restart diagnosis (dumps every 2e-5) | deterministic; one-cycle positivity failure at r = 9.1635 (inside WB zone, r_top 9.18), phi 61 deg, in a cold dense stream-fed layer (rho ~1e-2, v_phi 330) under hot thin ambient (rho 5e-7): contrast 2e4, shear 300 km/s over 3-5 cells; same class as RY Per t = 0.069 |
| gpu2/q11sm, q11_s45 | 12115722 | plaskett_env11 = env10 + fofc + nghost 3 + wb_rmax = R_acc | 0 FATAL to t = 0.045 (past 0.0353); dt 1.8e-6 -> 9.3e-7; FOFC firing almost all in floor/ambient cells r >= 9.34 |
| gpu2/q12sm, q12_s45 | 12115955 | plaskett_env12 (env10, wb_rmax = 8.52 = R_acc - 1.5 d_pen, nghost 2, no fofc) | 0 FATAL to t = 0.045 |
| gpu2/q12_a1 | 12115956 | env12 envelope only, 0.2 orbit | 0 FATAL to t = 0.0917 |
| gpu2/q12fsm, q12f_s45 | 12115957 | plaskett_env12f (env12 + fofc + nghost 3) | 0 FATAL to t = 0.045 |
| gpu2/q12f_a1 | 12115958 | env12f envelope only | 0 FATAL to t = 0.0917 |
| gpu2/q12_s1, q12f_s1 | 12116831 / 12116832 | 0.5 orbit (tlim 0.2291), 2 apu nodes | **cancelled while pending (user, 10-07)**; DeltaAI ran its copies |

Figures of q12_s45 / q12_a1 / q12f_a1 (maps, budgets) are in viper's `page/plaskett_test/` (not copied).

## Literature relied on (only what the notes cite; text dumps in viper lit/, list in lit_files.txt)
RY Per: Barai et al. 2004, Olson & Plavec 1997, Kolbas et al. 2014, Van Rensbergen & De Greve 2016,
Sudar et al. 2011, Etzel & Olson 1993, Peters & Polidan 2004, Richards & Albright 1999; stream scaling
Lubow & Shu 1975; Roche lobe Eggleton 1983. Candidates (CANDIDATES.md): Lechien et al. 2505.14780,
Schootemeijer et al. 2018, Wang et al. 2026, Vinciguerra+2020, Shepard et al. 2025, Harmanec & Scholz 1993,
and others listed there. Plaskett: **Wade et al. 2026 (arXiv:2609.25526)**, start-of-accretion quote
(M2 ~16 Msun, R ~9 Rsun, Mdot ~1e-4 Msun/yr; initial 18.2 + 16.8 Msun, P_i 3.69 d).

## Known pgen limitations
- Stage 1 (inner = surface): isothermal only; `thermo = adiabatic` FATAL there. The envelope mode is
  adiabatic (ideal gas) only.
- 2-D r-phi in a 4-cell reflecting theta band: no vertical stream structure; Mdot in code units is per
  band, so code Mdot scales with the phi width only (see NOTE section 2).
- The theta band shows 10-100x cell-to-cell structure at the q10 failure interface (PLASKETT.md); a true
  r-phi run needs nx2 = 1, which AthenaK does not allow with nx3 > 1.
- Cumulative flux integrals restart from 0 on a restart.
- scripts/ana_env.py shell radii are hard-coded for RY Per (R_acc 4.06).
- The comment block at the top of plaskett_env*.athinput is still the RY Per header (stale), and the
  env_amb_rho comment "(stream peak 1, rho_ph 20)" refers to the RY Per rho_ph.

## Open problems and the next steps viper would have run
1. q12_s1 / q12f_s1 to 0.5 orbit (now DeltaAI's), then the envelope-only drift check at 0.2 orbit and,
   if both arms are clean, choose wb_rmax = 8.52 with or without fofc.
2. Spin 7.2 twin with env_r_spin = R_acc (env7) once spin 1 runs 0.5 orbit.
3. Run long enough for a quasi-steady state (in t6/t7 ~90 % of the stream AM was still piling up outside
   R_acc at 0.5 orbit); report j_acc / j_K, Mdot_acc / Mdot_in over windows.
4. Resolution check in phi (1024/2048/4096) per DESIGN staged plan.
5. The stream density scale / width (NOTE section 2): fix the env_rho_ph derivation before any production.
6. Radiation (M1) on option 2, PLASKETT.md: ~1 week of setup (units, opacity, inner-flux BC, radiative HSE),
   cost ESTIMATE 5-9 h/orbit/node.
7. Core-change candidates, not done: per-cell WB switch under inflow/strong shear; positivity-preserving
   reconstruction limiter.
