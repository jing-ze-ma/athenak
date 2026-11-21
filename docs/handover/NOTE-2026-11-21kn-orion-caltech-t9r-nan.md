# NOTE 2026-11-21kn orion -> caltech: t9r NaN stage on the new settings — apply the kl fallback now, post logs

t9r (log.1) cycle 26418 t2s 2: 200 solves, FCONE/FWB momentum sums -nan, write-back closure inf -> BE stage fallback; cycle 26420 again
(neg 5955, closure 64). That is a NaN inside the solve with fwb+fcone at cfl 0.9.
1. t9r: per NOTE kl, fall back NOW to cfl 0.6 with the same keys (restart from the same backed-up rst if the run is not clean after
   26420; say which). If it recurs at cfl 0.6, go to flim (clean_stag_flim=true, fwb/fcone false) at cfl 0.6.
2. Post the t9r log from the restart to cycle 26425 (gz) to caltech_results/PROD1121b/t9r/k9_log.gz, and the restart rst name/cycle.
3. Post the first 50+ cycles of the other three runs (t9l62, m35mms, i1_m35mms_f1k) on the new settings as soon as you have them.
Do not change anything else.

-- orion (for the USER)
