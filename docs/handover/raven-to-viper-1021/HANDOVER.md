# HANDOVER 2026-10-21 ~20:00 CEST: Raven coordinator -> viper

The Raven coordinator session crashed some time after ~16:00. This note covers everything since then.
Full state: `raven_resume_1021.md` (the last Raven RESUME block, 05:00 to 16:00 entries) + worker NOTEs in this dir.
/raven/ptmp is not mounted on Viper, so every branch below is now PUSHED to the fork.

## Side branches (all opt-in, default bitwise, UNMERGED; pushed 10-21 ~20:00)
| branch | tip | what |
|---|---|---|
| gsdrag-1021 | 24810dad | `<rad_m1>/implicit_cs2_dacap = K` (+ realisable matched drag ref). BE+cs2+K2 = GO level (9.07k fast, 1.11x v11), gates PASS (see gsdrag_1021_NOTE.md) |
| pgrst-posfix-1021 | 2c53b9a2 | `implicit_hold` (pgen-held cells = Dirichlet rows in every implicit solve); design in pgrst_1021_NOTE.md. **NEVER BUILT/RUN**: the Raven build (interactive 31373149) was killed by the time limit at 17:52. Next: build, then arms h1-h4 (v12lzc hesdirk2 + cs2 + dacap + hold: stage vs BE steps, fast cells, budgets) + default bitwise gates. Held-energy face-flux booking is still TODO |
| pgrst-comb-1021 | a1c5d379 | k1_reset=pos + vimp_ewf(+Jacobian) + dacap + implicit_vimp_kill |
| pgrst-vimp-1021 | 33f3a46e | vimp_ewf, implicit_cs2_wbsolve (PARKED by USER) |
| pgrst-1021 | 2de67c53 | booking B / K1 (incorrect, do not merge) |
time2_posfix (e1612750) and time2_accept_fallback (5b525abb) are in the posfix history. Results: posfix does not cure cbg/cbh, and accept_fallback accepts 0 stages (S2 <= 0 is genuine).

## Production (no change without USER; fin51-r-fix = 077dda23)
- Raven, all pending: hgP145 31263626 (est 10-22 12:39), rsgP25 31316848 (14:16), acc3dP run12 31291344 (15:16, BE+matched v12lzc_be_gm), hgGlobR 31337226 (19:48).
- Viper: V2 12291385 / HG V 12278622 (he_apu 4de8168d_r71). Check their status.
- DeltaAI RSG E: CANCELLED by USER (NOTE 8c925750); 0-24.5 d kept on /projects/bivj. RSG E page scripts went to Caltech (d316b8a0).
- Caltech: H200 reproducibility fix = UCX_RNDV_SCHEME=put_ppln + GDR off (548404ec). RSG E fresh on 63a7d28b (bitwise = 4de8168d, -16.7 %/cycle) and HG global started (8dc031d7). Caltech accretor v12 chain CANCELLED (USER, 3bcec965).

## OPEN, not started by Raven
1. **TASK-2026-10-21-raven-accretor-vetgd.md (USER)**: test the vet_gd closure on the 3-D accretor (sp theta wedge nx2=58, walls at both theta faces). Check support -> short test vs vet_col (cost, dt, redos, budget, field diffs) -> template + md5 if GO. Nothing was run.
2. implicit_hold build + arms (above).
3. USER decisions pending: (a) accretor matched -> BE + conservative2 + dacap=2? (b) dacap in HG global (changes 1e-3 by design)? (c) DeltaAI dump cadence (moot now it is cancelled), Viper chains 6 vs 12 nodes, time2_k1_reset default.
4. TODO 10-22 (USER): tau-free central/half-range blend study (idort_f vs idort_a vs f vs H/J + Eddington-tensor anisotropy) on the IDORT paper tests + order grid. First item: thin waves (tau <= 0.1) show E order 0.66-1.6.
5. Pages (new account): IDORT page and others were republished. URLs are in Raven memory page-urls (IDORT Ve7rx, RSG g25 87Evb, Golden GajcP, RSG E DeltaAI LBASg, accretor D8q68, He giant THKEf). The Claude Docs page 'Shock-Supported RSG Atmospheres' still needs a manual export.

## Rules (USER)
Globals (RSG E / HG) stay on 4de8168d + hesdirk2 + conservative2; only the accretor is on 077dda23. Gates = current prod decks (RSG E clip N32, HG global clip N32, hg145, accretor v12lzc). One prod branch (fin51-r-fix). Every new feature needs >= 2 identical runs bitwise + continuous-vs-restart bitwise. Don't rerun known arms. Builds: incremental, on a CPU compute node.
