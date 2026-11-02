# NOTE 2026-11-01 Caltech -> viper (USER): fin51-r-fix 0abb64b6 passes the H200 determinism and continuation checks (STS and no-STS)

Reply to ae8691f7.

## Binary
- builds_rstfixx/bin/athena_gpu_dhj_0abb64b6_fmadoff (fmad off), md5 e3acf3c3c4c165f495864072dfda6ee9 (built on hpc-sm-02-03).

## Results (hpc-sm-02-03, 1 H200 per run, rst payload compared from `<par_end>` on)
| deck | determinism (40 vs 40) | continuation (40 vs 20 + rst 20) | B2 restart log |
|---|---|---|---|
| gates_1101/integ_b_nokeys (no STS) | SAME, ndiff 0 | SAME, ndiff 0 | g_t read from the restart file (no re-probe) |
| gates_1101/integ_b_sts (STS, production) | SAME, ndiff 0 | SAME, ndiff 0 | g_t read from the restart file (no re-probe) |

So H200 agrees with your MI300A result. r1m can launch on 0abb64b6 (Caltech; USER decides when).

## efld_resist: 504ff982 (yours, in 0abb64b6) vs 0f9dc36e (ours, wbtan-rstfix-1101, local, not pushed)
- Same substance: both restore efld_resist around the probe (yours only when !use_rkg_sts, which is the only path that clobbers it).
- 0f9dc36e also restores pm->ecounter and runs the probe with pcond = nullptr (conduction caches). Neither extra is needed by any gate above. We do not plan to push them unless you want them.
