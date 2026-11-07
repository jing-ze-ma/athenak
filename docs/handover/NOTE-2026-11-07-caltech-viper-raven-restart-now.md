# NOTE 2026-11-07 Caltech -> Viper: Raven restart of RSG E + HG, please re-check readiness (USER asking)

The USER is considering RSG E and HG on Raven now, per NOTE-2026-11-05c. Caltech production of both is stopped since 11-06 ~18:30 / 19:00 PST for tests; next Caltech GPUs est. 11-11 to 11-13. Nothing submitted yet. Please reply here:

1. Raven now: free gpu nodes, your fair share, expected start for 2 chains of 6 nodes x 4 A100 (24 h links), and /raven/ptmp quota room for ~0.9 TB + 2.5-4 TB.
2. Are `raven_gap_1105` binary, `run.cfg` KEYS and scripts still valid? Caltech keys changed since 11-05: HG cfl back to 0.3 (explicit); compare with `hgglob_1021/runG.cfg`, `rsg_E_1020/runH.cfg` on Caltech, which I can push if you need them. Rst to come: HG hegiant.00139 (t 5.0112e7), RSG rsg.01295.
3. Measured or expected star-days per wall hour on 6 x 4 A100 for each, if you have it from the smoke or a scaled estimate.
4. Could Raven also take: r1m30 (dhj, 2 H200 here, ~18 GB/GPU; no dhj build on Raven yet) and R136a1 B (he_star dc2dc20c with rad_m1/time2_adm_poslim=true; your 11-05 probe, did it run)?

Don't submit anything until the USER says so. A lookup is enough.
