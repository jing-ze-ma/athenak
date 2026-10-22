# TASK 2026-10-22 caltech -> viper (USER): make <rad_m1>/time2_accept_fallback default-on under hesdirk2; use it in accretor prod if all gates pass

Result (Caltech, RSG E restart rsg.00021 at 105 d, 1 node H200, 4 ranks, tlim 9.40e6 = 108.8 d; production onset = cycle 9260, 107 d):
 A d3df1401 harmonic (control)                      : 8 DT_RETRY after cycle 9200, dt -> 720-1065 s   (reproduces production exactly)
 C ff1f52bc harmonic + time2_vet_col=extrap2        : 10 after 9200, dt -> 458 s                      (extrap2 alone: no cure)
 E = C + time2_accept_fallback=true                 : 0 after 9200, dt = 1000 s (CFL) to 108.8 d       PASS
Cause of the 107 d retries: hesdirk2 stage 1 rejected with converged=1, min(E,T)=+1.5e-6 because the implicit_vimp / offdiag operator
positivity fallback fired in a Picard pass (fbk path) -> time2_reject = retry cuts dt. Your own comment (pgrst-posfix-1021): the BE redo
uses the same fallback operator, so rejecting only loses order. E and C both show an extrap2 restart blip (cycles 9073-9083, negative
E, BE redo) that recovers in ~50 cycles; E also had 4 isolated recovered retries before 9200.
Physics context: RSG photosphere / HG photosphere break into density clumps (dens/shellmed >100 or <0.01 in 1e5-3e6 cells at 1.0-1.08 R_ph,
all panels; HG saturated by 15 d) -- physical (RT-type), the solver must stay stable with them. The true RSG collapse (neg E, dt 2-5 s at
129.8 d) is being tested from rsg.00025 (E2 = E + time2_posfix=4; C2, D2 = ff1f52bc retry fix alone) -- results follow in a NOTE.

USER asks viper:
1. Flip time2_accept_fallback to default true when time_scheme = hesdirk2 (keep the key to turn it off). Gates as usual: bitwise for runs
   that never take the operator fallback (he RSG E det, accretor D_B1v6), C2==C and R==C with the key on, and the accretor / BSG A/B.
2. If ALL gates pass: use it in the accretor production (run12 / next link) and report.
Caltech productions RSG E + HG relaunch from the pre-onset rst with ff1f52bc (md5 cd998ac0) harmonic + extrap2 + accept_fallback=true
as soon as their allocations start (no tests in them). Still also waiting on the harmonic default flip (9c432f17).
