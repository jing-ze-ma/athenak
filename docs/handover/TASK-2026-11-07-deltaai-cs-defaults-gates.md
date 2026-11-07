# TASK 2026-11-07 viper -> DeltaAI (USER: "can you also do those on delta or deltaai"): cs default-candidate gates on hgb / rsb, ghx4-interactive

Thanks for 0a4be4102. Goal (USER): decide which opt-in keys become defaults for the RSG/HG productions, as fast as possible. Your interactive queue (~10 min) is the fastest place for these gates. Raven/viper run the same gates in their dev queues; first result wins, duplicates are fine. Budget: at most **2 interactive jobs x 1 h x 1 node (16 SU)** for this TASK.

## Inputs (USER copies from Raven to /work/nvme/bivj/jma20/cs_1107/)
- `/raven/ptmp/jinma/caltech_1107/rst/small/hgb_hegiant.00001.rst` 6,670,021,763 B, md5 596fe2d4525248df08f75f8bf84036aa
- `/raven/ptmp/jinma/caltech_1107/rst/small/rsb_rsg.00001.rst` 3,968,947,892 B, md5 4c1c93d040472e3b26e9f0e03e1f52e1
- `/raven/ptmp/jinma/caltech_1107/deps/raven_deps_1107.tar` 1,167,595,520 B, md5 426a2d8e106c022f3f97fbf36cb0c276 (tables; unpack, tree = Caltech $J)
- Table keys to repoint on every command line (absolute Caltech paths inside the rst): HG problem/he_ic_file, he_opac_table, he_planck_table -> hgglob_1021/files/hg145-prod-1018/{ic_hegiant_14p5M_mk.txt, rosseland_mesakap_afe2dd8da2bd.txt, planck_ext2_x0.0009496_z0.0215304_416a61196c64.txt}; RSG -> rsg_E_1020/files/rsg-E-1020/ic_rsg_g25_deep.txt, rsg_E_1020/files/rsg-g25-prod/{rosseland_mesakap_4c0e8b64d26a.txt, planck_ext2_x0.653848_z0.019757_d38e5526afef.txt}.

## Binary
Your `athena_hesnf_gpu_698b8241329e` (md5 136cf6db3ce6e79950c0e0ca7177c09a), plain `srun -n N --cpu-bind=cores` (the --gpu-bind=closest hang, 10-09).

## Keys
- BASE HG = Caltech production keys: `output2/dt=28800 output3/dt=28800 rad_m1/time2_vet_col=extrap2 rad_m1/time2_accept_fallback=true rad_m1/time2_posfix=4 rad_m1/implicit_work_fix=true rad_m1/cs_vertex_mono=true rad_m1/implicit_pos_floor_cons=true rad_m1/time2_auto_be=16 rad_m1/time2_k1_stiff=1.0 rad_m1/dbg_t2_admiss=20 rad_m1/vet_gd_twin=false rad_m1/vet_gd_symref=none rad_m1/vet_gd_twin_fuse=false rad_m1/vet_gd_twin_state=mean rad_m1/vet_source_noesrc=false time/dt_retry_host=false rad_m1/vet_gd_gl_nmu=0 rad_m1/vet_gd_geo=2 rad_m1/vet_rotate_avg=8 rad_m1/vet_rotate_advance=true mesh/cs_seam_blockend=extend mesh/cs_vertex_3cell=true mesh/cs_seam_rho_guard=0 time/cfl_number=0.3 `
- BASE RSG = `rad_m1/time2_vet_col=extrap2 rad_m1/time2_accept_fallback=true rad_m1/time2_posfix=4 time/dt_retry_t2_be=0 rad_m1/implicit_work_fix=true rad_m1/cs_vertex_mono=true rad_m1/implicit_pos_floor_cons=true rad_m1/time2_auto_be=16 rad_m1/vet_gd_twin=false rad_m1/vet_gd_symref=none rad_m1/vet_gd_twin_fuse=false rad_m1/vet_gd_twin_state=mean rad_m1/vet_source_noesrc=false time/dt_retry_host=false rad_m1/vet_rotate_advance=true mesh/cs_seam_blockend=extend hydro/merge_lat_blockend=exchange mesh/cs_vertex_3cell=true mesh/cs_seam_rho_guard=0 `
- r8-cs: `rad_m1/cs_cross_mono=true rad_m1/implicit_opac_newton=false mesh/cs_seam_rho_guard=4 rad_m1/implicit_predictor_tcap=1.1 rad_m1/implicit_fb_local=true rad_m1/implicit_fb_local_rounds=4 rad_m1/implicit_opac_lag_passes=3`
- T: `rad_m1/vet_gd_tangential=true`
- every run: outputN/bin_double=true for the bin outputs, hst on; rst output dt huge (rst only at the end).

## Job 1 (hgb), job 2 (rsb): ONE script each, all steps in sequence, 1 node, <= 1 h
1. R1: 1 cycle from the Caltech rst with BASE (writes an rst with the extended rad_m1 record; rst_state_ext is default on in prod-1107v; the Caltech rst lacks it).
2. G1 restart: from R1 with BASE+r8-cs+T: +4 cycles straight vs +2 -> rst -> +2; cmp the final bin dumps.
3. G2 layout: from R1 with BASE+r8-cs+T, 4 cycles on 1 GPU vs 2 GPUs (vs 4 if the blocks divide); cmp.
4. G3 determinism: from R1 with BASE, 4 cycles, twice, same layout; cmp.
5. If time remains (stop 5 min before the limit): short A/B from R1, same layout, equal cycles (as many as fit, e.g. 50-100): BASE / BASE+r8-cs / BASE+r8-cs+T / BASE+r8-cs+problem/rst_pert_erad=1e-10. Report per run: NOT POSITIVE, NOT ADMISSIBLE, BE redo, fallback/FBDIAG lines, m1_pos_flr floor counters (RSG: Caltech saw fb_local x2.6 floors), Picard mean, dt, s/cycle, GPU memory.
**Preflight first (USER rule after 5 jobs died at start tonight):** on the login node, check every command line: fresh start -> every key must exist in the deck; restart -> the <block> must exist in the rst input (absent keys are added). Script: `git show fork/bsg-files-1009:docs/handover/caltech-to-raven-1107/../../../` is not needed; the logic is 40 lines, copy it from viper's `/raven/ptmp/jinma/pool_1107/preflight.py` if the USER can bring it, else check by hand.

## Report (NOTE on this branch)
Per gate BITWISE / DIFF (which output, max rel diff), the A/B table, GPU memory per GPU, s/cycle, SU used. Stop rules: FATAL, NON-CONVERGED=[1-9] in the summary, NaN, rc != 0 -> stop and report.

## Not now
Productions on DeltaAI: not in this TASK (budget = ~one production for a few days; the USER decides). Holders: not now.
