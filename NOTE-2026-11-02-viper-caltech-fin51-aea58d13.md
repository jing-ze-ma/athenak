# NOTE 2026-11-02 viper -> Caltech (USER): fin51-r-fix = aea58d13b (follow-up: RT flux-hst restart cache + your probe extras)

fin51-r-fix aeafae97 -> aea58d13b = merge of followup-1102 (69e172aaa), same tree as 69e172aaa:
A. RT flux-history restart cache (rthst-1101): rt_Fb, rt_Qb, icut, mu0, fstar, albedo carried across a restart in three
   new implicit-ck restart slabs, so the dhj.user.hst RT flux columns (Lir_top Lsw_refl Lsw_abs Lrad_bot Etot_top
   Etot_bot Lsw_in) no longer read 0 in the first rows after a restart off the ck cadence.
   problem/flux_hst_rst default = flux_hst && ck_implicit (user); explicit true without ck_implicit warns, ignored.
   Restart files grow by 3 slabs (W121 C32: +26.7 MB). Old files load (gap rows stay 0 once); the old binary reads
   new files with the state bitwise. flux_hst_rst = false: rst bytes as before.
B. Port of the two extras of your 0f9dc36ec (credited): ProbeT saves/restores Mesh::ecounter and detaches pcond.
   Your efld_resist part not taken (fin51 already has 504ff9823, saved only when !use_rkg_sts). Please check that
   your wbtan-rstfix-1101 needs nothing else.
Effect: state bitwise unchanged; only the event log's first interval after a fresh start loses the probe's floor
hits (r1m bare: -11328 dfloor, -9596 efloor), and restart files gain the 3 slabs.

Gates: viper MI300A: 20+20 vs 40 rst data + mhd.hst IDENTICAL with STS and without; base 0abb64b6 vs new 40 cycles
mhd.hst IDENTICAL, with flux_hst_rst=false whole rst file IDENTICAL; W121 dhj.00641 +20 cycles mhd.hst identical.
Raven CPU: mhd *_cpu + rad dhj ck/srclim 34 passed; C8x112 150+150 vs 300 rst data + both hst IDENTICAL incl. RT
flux columns (first continued Lir_top 2.80e30, was 0). GPU user.hst Efloor/Mfloor/Efloor_rt differ 1e-16..5e-14
between ANY two GPU runs (reduction order), as before.

Please rebuild aea58d13b for r1m (or keep 0abb64b6 if a run is already going: physics bitwise the same).
