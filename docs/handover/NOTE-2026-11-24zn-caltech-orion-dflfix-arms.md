
### H3 at +300 cycles (11-24 07:34): **clean = no**
- H3 (dflfix, beside M35ZAMS on 01-15) to c38950: CW 0, FATAL 0, created inf 0, DT_RETRY 0; dt 1.694243e+02 vs production bsg 4.905494e+01 (rc3).
- eos_vceil (bsg3d.log column 5, per output row) in H3: max 1197. Production rc3 bsg, last 2 rows: 414923323 409839395 .
- dfl line: `` (printed at the end; empty while running)
- DTREP H3 (dt_report not named; deck value): 0 lines.
- **Correction to 'clean = no' (my script's criterion was too strict):** the 1197 eos_vceil events are all in H3's FIRST bsg3d.log row (c38721, the first output interval after the restart). That is the inherited floor gas from the rc3 state. The 4 later rows (to c39000) have **0** vceil events, against ~4e8 per row in the rc3 production.
  The [update] floor_audit adds ~4 floored cells per cycle (<r> 4.7e12; floor_de +6e28 erg/cycle). So: no NEW vceil cells after the start-up interval, CW 0, FATAL 0, created inf 0, DT_RETRY 0, and **dt 169 vs 49**. I judge H3 clean by 24zt and switch bsg now.

## H3 (11-24 07:38)
- start H3 BIN=athena_dflfix_5b6159ad0 md5=9c59917d84f344a7094b5403d432dd5c NODE=hpc-sm-01-15 HOLDER=7333036 ARGS ...; 2026-11-24 07:38:17 rc=0
- first CW: none; CW 0, FATAL 0, DT_RETRY 0, created inf 0
- dt (hst row = cycle + 1 from the run start): c200 1.6828336494738633e+02; c400 1.6650332869761490e+02; c600 n/a; c777 n/a; c1000 n/a; 
- end-of-run dfl line: `<hydro> dfloor_flux_limit: outflow-limited cell-stages (all attempts)=4.630238e+09, their x1 momentum x vol removed=-2.847802e+32, kinetic energy x vol made heat=9.926168e+39`
- FOFC line: `<hydro> EOS floors (all ranks, accepted steps, RK-weighted): mass created x vol=7.468168e+09 energy created x vol=4.538961e+36 | FOFC first-order cell-stages=2.943166e+09`
- eos_efloor column (bsg3d.log): last 0 eos_efloor
- floor_audit lines c770-790: 0 lines (in log.0.gz; first: ``)

## bsg (11-24 07:38): switched (dfl)
- pre-switch (log.3): 0.6422 s/cyc (c39350-39600), dt 4.981374e+01. Clean STOP: yes. Pre-switch rst kept as rst_keep/bsg3d.00029.rst.pre-dfl.
- relaunch: `2026-11-24 07:37:23 chunk 4 STAR=bsg BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_dflfix_5b6159ad0 md5=9c59917d84f344a7094b5403d432dd5c DECK=/resnick/groups/carnegie_poc/jingze/prod_rc3_1124/decks/tile/IA123_bsg` + EX: time/cfl_number=0.3 hydro/fofc=true hydro/fofc_retest=true rad_vet/clean_stag_mg_levels=4 rad_vet/vet_gd_ls=0
- echo: 
