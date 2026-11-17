# NOTE 2026-11-17 Orion -> viper (USER): push 18 CPU regression of b53e35565 vs push 17 5785ee51d

This answers TASK-2026-11-17-viper-orion-push18-regression.md.

**A: bitwise, with no movers.** Every reg run of b53e35565 with no keys is byte-identical to 5785ee51d: .bin data,
tab/hst and every CLEAN line. This covers the multi-rank decks, x1line, WB and the RT feed decks. No reg deck has a
run-time vlat_icut > 0 that 5eb7adf3f would move.
**B: no change from push 17.** maxrel is 0 on every deck, and 0 CLEAN WARNING / 0 DT_RETRY / 0 FATAL.
**Extrap arm: rc 0, 0 warnings, 0 unforced retries.** Conservation holds to round-off. The frac rule (1e4cfd388)
leaves rw bitwise. Separately, extrap on its own moves rw F1 by about 35 %, and that does not shrink with N (see the
table).
**Restart: bitwise.** That holds in all 9 A restarts and the 3 extrap restarts. The forced retry with extrap
(cycle 4, fac 1) is bitwise to the straight run on sp3, thin1 and hcyl.
**Timing (Raven, dedicated node):** flz_q +1 %, cs_closed -1 %, clean stage -1 % vs push 17. All are inside the
run-to-run spread.

## Binaries (Orion, git archive trees, gcc/13 + openmpi/4.1, Release, MPI+OpenMP, -ffp-contract=off, NMAX 4096, built_in_pgens)
| head | md5 |
|---|---|
| b53e35565 (push-18 candidate) | `851a7eb5e3c7c51573b03a23de4d595a` |
| 5785ee51d (push 17) | `941c8be8b938122fec3d2220521a4120` |

- Both were built in fresh dirs with kokkos d8e9af031.
- The runs used hold 231007 (A arms) and hold 231003 (extrap / rw arms), as `srun --overlap` steps with
  OMPI_MCA_hwloc_base_binding_policy=none.

## A: b53e35565 vs 5785ee51d (no keys)
- **All runs have rawdiff 0, bindata_diff 0 and CLEAN same.** That covers:
  - the cc ladder (hcyl, cbg, cb, shd, nltek1);
  - thin_p/r/f, dif, lat, pulse0, s88, s81 and rhse;
  - cs uni / closed / thin1 / thin1mb / thin1mbmr (2 ranks) / feed1;
  - flz_q and flz_Q;
  - rw x4;
  - RT hcyl, sp3mb, sp3mbmr (2 ranks), feed, feedx1 (x1line), feedmr (2 ranks), thin1, rot1s and csCnr.
- WB fv1 / fv / lld / cart: rh.user.hst are byte-identical.
- No mover, so there is no first differing file to report.

## B table (b53e35565; every value equals push 17 5785ee51d)
| deck | maxrel vs P17 | CLEAN lines | metric |
|---|---|---|---|
| cc hcyl cbg cb shd nltek1 | 0 | same | scores 1.2061e-02 1.1003 0.45999 1.1805e-02 9.0716e-05 |
| thin_p / thin_r / thin_f | 0 | same | L1(E) 1.289e-01 / 1.179e-01 / 2.942e-02 |
| dif / lat | 0 | same | L1(E) 1.536e-02 / 2.915e-01 |
| pulse0 | 0 | same | dE +8.882e-16 |
| s88 / s81 | 0 | same | steady pass 1: 2.578e-07 |
| rhse | 0 (hst) | same | vcs_max 1.85058e-01, it 806 |
| cs_uni / cs_closed | 0 | same | uni max\|E-1\| 0; closed sum E V rel 0.00e+00 |
| cs_thin1 / thin1mb / thin1mbmr | 0 | same | thin1 L1(E) 3.891e-02 |
| cs_feed1 | 0 | same | dE rel 0.00e+00, dLz -5.84e-17 |
| flz_q | 0 | same | dLz 5.87e-16, dE 8.40e-14, it 7909 |
| flz_Q | 0 | same | dLz 1.25e-15, it 711 |
| rw x4 | 0 | same | |
| WB fv1 / fv / lld / cart | hst identical | | vcs_max 1.21364e-09 / 2.96263e-09 / 4.75321e-09 / 7.38032e-09 |

- 0 CLEAN WARNING, 0 DT_RETRY and 0 FATAL in every A run of both builds.

## Extrap arm (`rad_m1/clean_fv_lag = extrap` written into the deck copies; b53e35565)
| deck | rc / W / DT / F | maxrel vs no-extrap | "clean_fv_lag = extrap:" one-pass / swept on (largest miss/change) | CLEAN its (vs no-extrap) | dE / dP |
|---|---|---|---|---|---|
| sp3 (sp fv1, force_reference = clean, 16 cyc) | 0/0/0/0 | 3.7e-5 | 14 / 0 (0.37) | 275 (270) | |
| flz_q (sp fv1, 40 cyc) | 0/0/0/0 | 6.9e-1 (chaotic deck, as with yres in push 17) | 5 / 49 (3.2) | 8959 (7909) | dLz 5.87e-16, dE 8.29e-14 (8.40e-14) |
| cs thin1 (RT deck, 16 cyc) | 0/0/0/0 | 2.7e-12 | 14 / 0 (3.2) | 68 (68) | |
| cs thin1 (reg deck, 100 cyc) | 0/0/0/0 | 7.8e-16 | 98 / 0 (**3.2e+273**) | 118 (119) | L1(E) 3.891e-02 (same) |
| cs feed1 | 0/0/0/0 | 4.3e-3 (hydro_u) | 7 / 1 (1.7) | 260 (257) | dE rel 2.22e-16, dLz -5.84e-17 |
| hcyl (Cartesian IDORT beam, clean_rays = fv1) | 0/0/0/0 | 6.4e-16 | 14 / 0 (**3.5e+10**) | 267 (267) | |
| rw N32 +-1e-3 (fv1 + vet_gd keys, see below) | 0/0/0/0 | E 2.5e-6, **F1 3.6e-1** (odd 3.6e-1) | 42 / 21 (1.1) | 497 (363) | |
| rw N64 +-1e-3 | 0/0/0/0 | E 5.0e-6, **F1 3.4e-1** (odd 3.5e-1) | 78 / 49 (1.4) | 861 (531) | |

- The forced arms (below) each show exactly one `DT_RETRY ... reason=forced`.
- **rw: extrap's change does not shrink with N.** F1 changes by 36 % at N32 and 34 % at N64, and E by 2.5e-6 and
  5.0e-6.
  - With `clean_fv_lag_frac = 0.0` (the push-17 rule only), rw is bitwise to the default frac 0.5. So this is the
    clean_fv_lag_tol rule (1e-4 of the largest intensity), not the new frac rule.
  - My reading: on a 1e-3-amplitude wave, 1e-4 of max I is about 10 % of the perturbation, so the tolerance is loose
    for linear-wave decks.
  - The no-extrap baseline is itself fv1 + the vet_gd ray set copied from hcyl. Against the default vet_sc rw (A18),
    that baseline moves F1 by 4.7e-2 (N32) and 3.5e-2 (N64).
- **Cosmetic: the summary's "largest miss/change" blows up when the inflow barely changes.** It reads 3.2e+273 on cs
  thin1 and 3.5e+10 on hcyl, where the miss/change denominator is close to 0.
  - Those builds were accepted by the tol rule anyway.
  - A floor on the change (for example tol x max I), or reporting the ratio only for builds the tol rule did not
    accept, would make the ratio meaningful.
- **Pre-existing, not new in push 18: `clean_rays = fv1` on the rw deck segfaults in VgdRotate.**
  - The rw deck is Cartesian vet_sc with no vet_gd_* keys. The call chain is CleanSolve -> CleanClosure -> CleanRays
    -> CleanFvSets -> VgdRotate.
  - The segfault is the same on 5785ee51d, at cycle 0.
  - Adding hcyl's vet_gd_* keys makes it run.
  - Under the "unsupported combinations must FATAL" rule this should become a FATAL with a clear message.

## Restart
- **A:** 9/9 RT restarts give `hydro_u same m1 same` on both builds. These are hcyl, sp3mb, sp3mbmr, feed, feedx1,
  feedmr, thin1, csCnr and rot1s (RTN).
- **Extrap:** 16 cycles vs 8 + restart + 8 is bitwise on sp3, thin1 and hcyl. The hcyl fv1 base without extrap is
  bitwise as well.
- **Forced retry + extrap:** the deck copies get `dt_retry = true`, `dt_retry_views = clean`,
  `dt_retry_force_cycle = 4` and `dt_retry_force_fac = 1`.
  - Their 16-cycle runs are bitwise to the straight extrap runs on sp3, thin1 and hcyl (hydro_u same, m1 same).
  - Their own restarts are bitwise too.
  - With the forced retry the summary counts one more build: 15 vs 14 one-pass.

## Timing (Raven job 32455624, one dedicated node ravc3029, Xeon 8360Y, `general --exclusive`)
- **Setup:**
  - Defaults, no keys.
  - flz_q: 40 cycles, OMP 8. cs_closed: OMP 6.
  - 3 rounds in the order a: P17, P18; b: P18, P17; c: P17, P18.
- **Binaries:** built in the job from the same git archives; md5 fe226669 (5785ee51d) and 1233e040 (b53e35565).

| build | flz_q s/cycle | cs_closed s/cycle | cs_closed clean stage steps (s) |
|---|---|---|---|
| push 17 5785ee51d | 2.659 (2.694 2.623 2.660) | 2.203 (2.185 2.214 2.210) | 3.416 (3.388 3.405 3.457) |
| b53e35565 | 2.676 (2.680 2.681 2.668) **+0.6 %** | 2.191 (2.202 2.137 2.234) **-0.5 %** | 3.396 (3.412 3.299 3.479) **-0.6 %** |

- All changes are inside the ±1-2 % spread. CleanFvRows shows no CPU cost or gain on these decks.
- The timing runs of the two builds give identical .bin md5s on both decks.

## Files (Orion: /orion/ptmp/jinma/clean_1115/push18_reg)
- Source trees and builds:
  - src_p17, src_p18, build_p17, build_p18;
  - athena_5785ee51d and athena_b53e35565 (copies of the binaries).
- Harness and arms:
  - reg/ is the push-17 reg harness copy.
  - The arms are A17, A18, X18 (extrap), R18 (forced retry), B18 (hcyl fv1 base), F18 / Z18 (rw fv1, rw extrap
    with frac 0).
  - The scripts are runA.sh, xarm.sh and tools/{cmpx.py, pd.py, rwcmp.py}.
- Raven timing: /raven/ptmp/jinma/stag2/p18tim (job18.sbatch, node18.log). Its runs are in
  /raven/ptmp/jinma/stag2/p16bis/reg (*BP18*).
- All .rst dumps were deleted (Orion and Raven).

-- Orion
