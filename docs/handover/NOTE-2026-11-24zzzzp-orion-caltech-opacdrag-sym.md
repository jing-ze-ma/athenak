# NOTE 2026-11-24zzzzp orion -> Caltech: opac_drag pairs (USER: ON) and the sweep question (sym)

T2 binary (7d7c4a71), idle nodes only (H0F / H1F of zzzzo first, then U4C, then these). All arms: lat_sums=true tie_gid=true thinf=false, vet_gd_allow_clamp=false, vceil=3e8, mg 4 for m35mms / 2 else, cfl 0.3, as in the U arms.
1. **PD (opac_drag):** per star UA' = `rad_vet/clean_stag_opac_drag=false` vs UD = `=true`, nothing else different; same rsts and cycles as U1-U4 (bsg, M35ZAMS, m35mms, t9l62 outburst). U1C / U2C / U3C / U4C are UA' where they were run with opac_drag unset (old default false): say so and reuse them instead of rerunning. Report and harmless rules as 24zzzzl (cost <= 1.05, L_top <= 0.5 %, no extra warnings / BE / retry).
2. **PS (sweep, solver only):** m35mms (mg 4) and t9l62 outburst, 100 cycles: `rad_vet/clean_pc_sweep=rb` (production) vs `=sym`; UC settings otherwise. Report Krylov its per stage, its/solve, passes, s/cyc, dt, CW / BE / retry, hst difference (expected solver-tolerance level). On the cs stars (HG / RSG fresh, old binary) sym was neutral: only these two sp stars have not been measured.
Do not run anything else on these keys.

-- orion (coordinator)
