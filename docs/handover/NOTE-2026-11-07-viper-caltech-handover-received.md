# NOTE 2026-11-07 viper -> Caltech: handover received; what runs where; please read before any Caltech restart

Re: HANDOVER-2026-11-07-caltech-to-raven.md, TASK-2026-11-07-bsg-novstar-raven.md (0caa4f54a).

## 1. USER decisions this evening (11-07, viper session)
- "use the latest updated branch for all runs" -> **all Raven/viper runs use `prod-1107v` 698b82413** (pushed to fork) = prod-1107 ee2b817d4 + vstarfix-1106 3b1d3dc19, clean textual merge (git merge-tree). The vstarfix keys (he_bc_linhao_fis, he_ic_linhao_seed_group, doppler_vel = vstar_nbe) are read only when named. This overrides your sec 1.3 recommendation (exact old commit for the first link): compare_hst vs your hst is now a tolerance check, not ~1e-12.
- "queue runs on both machines and see which start first": every run is queued on Raven (A100) AND viper (MI300A); first to start claims it (atomic mkdir under /raven/ptmp/jinma/bsgnv_1107/claims and defaults_1107/*/claims), the twin exits.
- "prioritize bsg", then "make those tests first. the goal is to figure out what to make default that can be used for prod runs. make those tests as fast as possible" (r7 / full tensor; sec 4).

## 2. Build finding you need
**Plain prod-1107 (ee2b817d4) does NOT compile with CUDA 12.6** (Raven): he_star_m1.cpp `ProblemGenerator::UserProblem` lambda #12 (4-D par_for, ~28 captures) -> "nvcc internal error: unexpected number of captures in __host__ __device__ lambda". prod-1107v compiles (its he_star_m1.cpp equals vstarfix's apart from one call line outside the lambda). If you build prod-1107 on CUDA 12.9 and it works there, it is a 12.6 nvcc limit; otherwise reduce the captures. Log: /raven/ptmp/jinma/builds/logs/he_a100_nofma_b_ee2b817d_1107_191533.log.
Binaries (prod-1107v): Raven he `athena_he_a100_nofma_b_698b8241_fmadoff` md5 ba86b8db65905bb67ab5daf1cbae2b38 (CUDA 12.6, AMPERE80, -ffp-contract=off --fmad=false); viper he `athena_he_gpu72_698b8241_p1107v` md5 79361078fc9f580d53479673b4235672 (ROCm 7.2); Raven dhj `dhj_a100_nofma` 698b8241 (W121). BSG reference = your prod_lhxo code: `athena_he_a100_nofma_3b1d3dc1_fmadoff` md5 ab8583d94d0af632a0a5ad98e677b664.

## 3. BSG no-vstar (your TASK) - queued, prod-1107v, deck = your novstar deck (diff vs bsgLI_lhx_exact_o = exactly the 2 intended lines; 4 tables md5-checked)
/raven/ptmp/jinma/bsgnv_1107 (viper twin /viper/ptmp2/jinma/bsgnv_1107). 1 node, 4 ranks x 4 blocks.
- STEP 0 gates (Raven gpudev): restart 2+2 vs 4 (tensor on), layout 8x2 vs 4x4, and **code switch: 4 cycles prod-1107v vs 3b1d3dc19 must be bitwise** (bins; rst differs by design: rst_state_ext default on).
- Arms from t = 0 to 2 d: A (deck: vstar off, tensor on), B (tensor off), Bn (B + he_ic_linhao_amp 0.05 x (1+1e-10), noise), **C = A + sp-generic r8 keys** (implicit_opac_newton=false, implicit_predictor_tcap=1.1, implicit_fb_local=true rounds 4, implicit_opac_lag_passes=3).
- Long run 0 -> 57.4 d, tensor on, queued NOW (2 x 24 h links) so it holds its queue place; if STEP 0 fails it is changed in place to vet_gd_tangential=false before it starts (your fallback). Energy gate during the run as in your TASK.
- Raven estimated start ~5 h from now (was ~25 h at 19:20 CET).

## 4. Defaults tests (sec 4) - running / queued
- R136a1 sp (rst 8, both machines): restart + layout bitwise with r8-sp + tensor; A/B control / r8-sp / r8-sp+tensor / noise (problem/rst_pert_erad 1e-10); dt collapses, posfix, NOT ADMISSIBLE, BE redo, energy budget. rad_m1/time2_adm_poslim=true named everywhere.
- He 128^2 sp wedge energy test with vet_gd_tangential on the compact face form (your T2a), viper.
- HG/RSG cs (r8 incl. cs_cross_mono + cs_seam_rho_guard 4, + tensor): restart/layout gates and A/B on hgb/rsb as soon as the small states arrive (transfer order: RSG, HG, then small; we asked the USER to send hgb/rsb in parallel).
- Proposal for defaults (USER decides after results): r8 without time2_k1_redist (your advice kept), tensor only if the He sp energy test passes. RSG: fb_local floors x2.6 still unexplained -> RSG without fb_local unless explained.

## 5. RSG / HG / R136a1 / r1m30 production
- NOT submitted. USER wants r8 + full tensor on RSG/HG -> waiting for the sec 4 gates. run.cfg on Raven updated (prod-1107v, table paths repointed, HG time/cfl_number=0.3 added).
- Node counts: Raven A100 40 GB -> RSG and HG each 6 nodes x 4 (1 block/GPU; memory-bound). Viper MI300A 128 GB -> RSG 4 nodes (8 x 3 blocks, your layout), HG 6 nodes (12 x 2, your layout). A 3-node RSG memory test (12 A100 x 2 blocks) is auto-queued when rsg.01295 arrives.
- r1m30: not queued (USER WB decision pending, your sec 2.4 (a)/(b)).

## 6. Please (Caltech / USER)
- Keep prodG/STOP and prodH/STOP.
- Before R136a1 continues anywhere off Caltech: create r136a1_1104/runB/prod/STOP or scontrol hold 6264757 (r136B2 would fork the history).
- If you push more to prod-1107 or the branches, say so in a NOTE: we build from prod-1107v and will re-merge.

## 7. FASTER TRANSFER (USER: "also about faster transfer")
Measured on Raven: the running rsync writes a flat **9 MB/s** (rsg.01295 partial 16.4 GB at 19:47 CET). Flat ~9 MB/s on a transatlantic RTT is the per-channel ssh window limit, not the link. At that rate the rest (RSG ~29 GB, HG 97.5 GB, small states ~30 GB) takes ~5 h.
Fix: several byte ranges of one file in parallel, each its own ssh channel. Over the existing ControlMaster this needs **no new password/OTP**. Expected ~3-5x with 8 streams (unmeasured; if it does not scale, open 2-4 extra masters, see below).

`pcp.sh` (run on the Caltech login node inside tmux, from $J):
```
#!/bin/bash
# pcp.sh FILE RAVEN_DEST [N=8]: N parallel byte ranges of one file over the open ssh master; writes RAVEN_DEST in place
F=$1; D=$2; N=${3:-8}; SZ=$(stat -c %s "$F"); MB=$(( SZ/1048576/N + 1 ))
ssh raven "mkdir -p \$(dirname '$D')"
for i in $(seq 0 $((N-1))); do
  dd if="$F" bs=1M skip=$((i*MB)) count=$MB iflag=fullblock status=none \
    | ssh -o Compression=no raven "dd of='$D' bs=1M seek=$((i*MB)) conv=notrunc iflag=fullblock status=none" &
done; wait
echo "local $SZ remote $(ssh raven "stat -c %s '$D'")"
```
Order (priority = the defaults tests, then production):
```
R=/raven/ptmp/jinma
./pcp.sh prod_1103b/he/runs/hgb/rst/hegiant.00001.rst $R/caltech_1107/rst/small/hgb_hegiant.00001.rst.part 8
./pcp.sh prod_1103b/he/runs/rsb/rst/rsg.00001.rst      $R/caltech_1107/rst/small/rsb_rsg.00001.rst.part 8
./pcp.sh hgglob_1021/prodG/rst/hegiant.00139.rst        $R/raven_gap_1105/hg/rst/hegiant.00139.rst.part 8
```
- Write to `.part` names; viper checks size + md5 (compute job) and renames. Sizes/md5s are the ones in HANDOVER sec 0.
- The running rsync chain: let it finish RSG (rsg.01295 is ~half done), then Ctrl-C it before it starts HG (otherwise it re-sends HG). R136a1 rst 5 and BSG rst 58 are low priority (R136a1 tests use rst 8, already here).
- If 8 streams give < 2x: extra TCP connections = extra masters, one OTP each:
  `ssh -o ControlPath=~/.ssh/cm2-%r@%h:%p -o ControlMaster=yes -o ControlPersist=12h -fN raven` (and cm3, cm4), then give each stream `-o ControlPath=~/.ssh/cmK-%r@%h:%p` round-robin.
- If both sides offer Globus endpoints (Caltech HPC and MPCDF): that would be faster still and needs no open terminal; we have not checked whether MPCDF has one for /raven/ptmp.
