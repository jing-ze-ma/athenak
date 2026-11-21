# NOTE 2026-11-21gx orion -> caltech: FOUND — CUDA 12.2 miscompile. Validation on B_v203n (CUDA 12.9), 4 arms in parallel

Verdict: B_v203n (e688cb7bc, CUDA 12.9.0 + g++ 13.2.0, md5 8d3f469f) on H200 STREAMS, bitwise the CPU/A100/GH200 numbers
(pass 1 4.020e-06, NOT converged after 50, F/c 1.3589e+42). B_v203 (CUDA 12.2.1 + g++ 11.5) traps, with arena off,
initcheck 0 and racecheck 0 -> deterministic codegen of the CUDA 12.2.1 toolchain. From now on: NO CUDA 12.2 builds;
every Caltech binary is built with CUDA 12.9 + g++ 13 (ccache must not reuse 12.2 objects: separate cache dir or disabled).
All B_* binaries built with 12.2 are retired.

Please run in parallel on 4 free nodes (former production nodes OK; one run per node, test nodes otherwise), all B_v203n:
1. BV_bsg: the UNCHANGED production deck bsg_wngl_fresh (d4ab303f), fresh, tlim 1.12e6 (day 13); hst, raycut, m1prof +
   m1_vet at 0.5 d, 6 d, 10 d and the end. No rst except a final one, deleted after.
2. BV_bsgtau: bsg_diag_tau1000 (aaed058b), same settings.
3. BV_t19: the unchanged t19l64_wngl_fresh (22d054ad), fresh, tlim 2.0e4 s; m1prof at 3 h and the end.
4. BV_cs: cubed-sphere check — the production RSG and HG restarts (as V203L_rsg / V203L_hg, v20.3 keys), 4 cycles each
   with B_v203 AND B_v203n (same node sequentially): report user.hst L_top/L_in and hydro KE rows of both, and the
   clean_init/clean summary lines, so we know whether the old HG/RSG H200 results were affected.
Also finish the B_v203r (12.2, no ccache) tiny smoke as the ccache control. Results to caltech_results/BV/.
No production launch (that needs the USER in your session).

-- Orion coordinator
