# NOTE 2026-10-19 caltech -> raven: 3-D accretor FRESH run with tracers (reply to TASK-2026-10-19-caltech-acc3d-fresh-tracers)

USER asked for this run. Status: SUBMITTED (build -> check -> 6 links, afterok chain). Nothing has run yet.

- Code: origin/acctrc-1019b @ 4960ac49, PROBLEM=ry_per_accretor, nofma recipe as v11 (-ffp-contract=off, HOPPER90, compiled on a P100 node).
  Binary: /resnick/groups/carnegie_poc/jingze/builds_accr/bin/athena_gpu_ry_per_accretor_4960ac49_nofma (md5 in build.5059464.out; follows in a later NOTE).
- Deck: docs/handover/acc3d-prod-1018/acc3d_prod_v11amwtr.athinput.in (git copy md5 5f372597, has @FILES@; Raven original 9fca5608 differs only by the paths),
  @FILES@ -> /resnick/groups/carnegie_poc/jingze/acc3d_trc_1019/files -> acc3d_prod_trc.athinput. Column / Rosseland / Planck md5 5b813886 / d81b7a46 / 99d668f6 verified.
  Layout 2 H200 (1 node), 8 MeshBlocks (4 per GPU), as v11.
- Check (job 5059465, 2 H200, afterok build): C, C2 static to t 0.002; R = C restarted from the t 0.001 rst; hst every cycle (dcycle=1 in a check-deck copy,
  acc3d_check.athinput) incl. ryper.trc.hst, final bins compared without headers. Exit 1 on any DIFF in C2-vs-C, or R-vs-C hydro.hst / trc.hst / bins.
  R-vs-C ryper.user.hst diffs are TOLERATED and reported per column with the offset range (cumulative columns, constant offsets after restart, Raven 39e9b9a6;
  the v11 check 5008035 failed on exactly this, which left v11 links DependencyNeverSatisfied).
- Production links (afterok chain on the check, 12 h each, WT 41999 s): 5059466 5059467 5059468 5059469 5059470 5059471.
  acc3d_chain.sh unchanged (md5 1ad77aed = git copy): FRESH from t=0, static phase (problem/stream=false) to t 0.05 in static/, then the deck as is in stream/.
  Stop on FATAL / Picard DIVERGED / summary NON-CONVERGED=[1-9] (not the per-pass "NON-CONVERGED after" lines) / NaN / rc != 0.
- Start estimate: none yet (build 5059464 pending on Priority, Slurm gives no START_TIME; gpu partition 111 running / 171 pending).
- Old v11 links 5007430-35 (held, DependencyNeverSatisfied) are left for the user to cancel; acc3d_v11_1019 untouched.

Paths (Caltech): run dir /resnick/groups/carnegie_poc/jingze/acc3d_trc_1019 (JOBS.txt, chain.txt, starts.txt, check/, static/, stream/);
hst stream/ryper.*.hst (user, hydro, trc); hydro_w dumps stream/bin/ryper.hydro_w.NNNNN.bin every 0.005 (every 4th = 0.02).
