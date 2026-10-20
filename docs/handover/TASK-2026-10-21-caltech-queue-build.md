# TASK 2026-10-21 raven -> caltech: queue estimates now + pre-build the fix tip (USER: move tests to Caltech if faster)

1. Reply within this poll (NOTE caltech->raven) with the CURRENT start estimates (`sbatch --test-only`) for:
   (a) 1 H200, 30 min; (b) 4 H200 / 1 node, 30 min; (c) debug QoS 30 min (CPU); (d) 1 CPU node 2 h (expansion).
2. Build now (no runs): branch `t2befb-ncfix-1cb` @ 4de8168d (= train tip 1cb24f82 + the cbg/cbh BE-fallback fix; not yet the
   prod tip; fin51-r-fix is pushed after Raven's gates). Same recipes as your 3896d529 builds: RSG he H200 bin and the
   accretor nofma H200 bin (pgen ry_per_accretor). Report md5 + `readelf -p .comment`.
3. Do NOT start or change any production (RSG E 5126654 stays as is; accretor stays cancelled until the v12lzc NOTE).
If (1a)/(1b) start in < 20 min, Raven will send short GPU test TASKs (accretor v12lzc screening arm, gate arms) to race Raven.
