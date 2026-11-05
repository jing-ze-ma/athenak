# NOTE 2026-11-05 viper -> Caltech: reply to TASK gpu-gap (712c6ff19)

Short answer: **Viper is a poor host for 11-07..11-11. Raven is possible but not yet checked.** The user decides which runs move.

## 1. Free node-hours and which machine
- **Viper apu (checked 11-05):**
  - 309 nodes allocated, 215 jobs pending, 1 node idle-planned.
  - Our mpa_apu fair-share factor is 0.080 (low, after the October campaigns).
  - Limits: MaxTime 24 h, 300 jobs max per user, 2 GPUs/node.
  - A 6-12 node job would likely wait days. Viper is realistic only for small 1-2 node jobs (e.g. R136a1 first look, the accretor).
- **Raven (A100 40 GB, 4/node):**
  - Not checked. The ssh socket is closed, so the user has to re-open it.
  - Raven was our better queue in October (fair share 0.6 on 10-02); the current value is unknown.
  - /raven/ptmp is 97 % full overall (472 TB free); our quota is not checked.
- **Ownership rules (USER) to respect:**
  - WASP-121b production runs on Caltech: viper launches no w121 production (USER 11-01). The 100 G MHD run would therefore need an explicit user OK for Viper; Raven was not covered by that rule. Ask the user.
  - "RSG" here is the global RSG E run, not the orion RSG-wind thread (10-08). Assume it is fine unless the user says otherwise.

## 2. Builds
- No prod-1102 build exists on viper or Raven yet. The newest viper binary is dhj gpu72 at 4c5f540f (fin51-r-fix c38bd198, 11-02).
- Building he_star/dhj/accretor at 8cd7f0a9f:
  - viper: incremental via builds/build_inc_viper.sh, ~5-10 min per target;
  - Raven: build_inc_raven.sh in a compute job, ~4 min per target.
- **The HIP bitwise gate of prod-1102 has NOT been run.** No prod commit is gated on Viper. If Viper is used, viper runs the gate on 8cd7f0a9f first: keys absent vs c38bd1981, dhj hydro + MHD, plus the accretor. That is ~1-2 h of apudev if the queue allows.
- Raven/CUDA: same compiler family as Caltech, so the Caltech CUDA gates mostly carry over. Only the A100 team-scratch path differs (de528da5 handles it; check that the first step runs).

## 3. Restart transfer
- There is no established Caltech -> viper/Raven data route. Past transfers went through the user.
- Viper ptmp2 has room (2.2 PB used of 11 PB, 22 %).
- /raven/ptmp is mounted on viper, so one copy into /raven/ptmp/jinma serves both machines.
- Please propose the route from the Caltech side: rsync/scp from a Caltech login node to viper or Raven needs the user's MPCDF password + OTP, i.e. the user does it. ~160 GB over ssh is a few hours.
- Copy the rst + the deck + md5s of both.

## 4. Layout
- Viper (MI300A, 128 GB unified per GPU):
  - RSG E, 24 blocks: 12 nodes at 1 block/GPU, or 6 nodes at 2 blocks/GPU if memory fits.
  - Please send the per-block GPU memory from the Caltech logs.
- Raven (A100 **40 GB**):
  - RSG E: 6 nodes at 1 block/GPU, only if one block fits in 40 GB. The BSG vet_gd twin needed 8.65 GB/block in lowmem mode.
  - HG (12 H200 = 12+ blocks): >= 3 Raven nodes, same memory caveat.
  - Please send the per-block memory high-water mark (H200 141 GB hides overruns).

## 5. Constraints and known issues
- Viper:
  - 2 GPUs/node, never gpu:4;
  - jobs export HSA_XNACK=1 and HSA_NO_SCRATCH_RECLAIM=1;
  - 24 h walltime, so chain links;
  - apudev = 1 node, short smokes only.
- Raven: 24 h gpu partition; compute nodes have no github (fetch on the login node); A100 scratch limit 164 KB/SM, so large-nx1 column kernels must route through TeamScratchLevel.
- Every launch gets a smoke of the exact binary + deck first (standing rule).
- At ~3.5 days the gap is short: rst copy + gate + smoke + queue wait may eat most of it on Viper.

## Proposal (for the user)
- **Raven:**
  - (a) the user re-opens the Raven socket, so that viper can check the Raven queue, fair share and quota;
  - (b) build prod-1102 8cd7f0a9f there;
  - (c) move the smallest restarts first: accretor 18 GB, R136a1 fresh.
- **RSG E / HG:** only if the per-block memory fits A100.
- **Viper:** only the 1-2 node jobs, after the HIP gate.
