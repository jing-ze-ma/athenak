# TASK raven -> DeltaAI (2026-10-18): He giant 14p5M -- baseline twin of Test B + background analysis of both (tests only)

Thanks for 93246a60. Your KE_int flag (x7 between 2 and 3 d in Test B) needs a baseline to decide floor effect vs convection onset.

1. Test B0 (1 GH200, ONE link, ~5 min): deck `docs/handover/hg145-prod-1018/hg145_mk_nobh.athinput.in` (baseline floor/sponge:
   dfloor 1e-18, pfloor 1e-11, he_sponge_dmax 2e-14; IC ic_hegiant_14p5M_mk.txt), SAME binary 4c494982 and EXACTLY the Test B keys:
   `mesh/nx1=704 meshblock/nx1=704 mesh/nx2=16 mesh/nx3=16 meshblock/nx2=16 meshblock/nx3=16 problem/he_seed=0.0
   problem/mlt_closure=adaptive time/tlim=2.592e5 output1/dt=300 output2/dt=4320 output3/dt=4320 output5/dt=1.728e5 output6/dt=1e30`
2. Background analysis on B and B0 with `bg145.py` (same dir; needs vis/python/bin_convert.py on PYTHONPATH via VISPY=<repo>/vis/python):
   `BASE=hegiant LOGF="<your run log glob>" python bg145.py 7.95704e12 2.22652e13 <rundir> 4`  (R_ph and sponge r0 in cm).
   Post the full stdout of both (it is a table per dump: per radial band <v_r> km/s and max|rho/rho0-1|).
3. For both: KE_int vs time (hst) at 0.25-d steps, and which radial band carries the KE growth after 2 d
   (shell-mean v_r and rms |v| per band from the last 3 dumps: r/R_ph 0-0.2, 0.2-0.5, 0.5-0.9, 0.9-0.99, 0.99-1.01, 1.01-1.03, 1.03-1.1, 1.1-1.5).
Report: NOTE-2026-10-18-deltaai-hg145-base.md. Nothing else to run.
