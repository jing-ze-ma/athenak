# TASK for Caltech: WASP-121b 3x-solar hydro production, fresh start to rot 300 (duplicate of viper w3x)

**User 09-30:** run the 3x arm on Caltech too; viper keeps its copy queued (job 12030356, est. start ~07:35 CEST 09-30).
**Before submitting, report the Caltech 1-node start estimate** (`sbatch --test-only`) in a NOTE; whichever machine
starts first runs to the end, and the user cancels the other copy (tell the user in the NOTE when yours starts).

## What 3x is
Third metallicity arm next to 1x and 10x (both at rot 300 on viper and Caltech). Observed day-night transport sits
between 1x and 10x (viper synth_rot300 and synth_rot300_10x: eps 0.35 / 0.13 vs obs 0.25). See
docs/handover/NOTE-2026-09-30-w121-3x.md (setup, grid nx1 76, IC, smoke) for everything.
- ck/CE tables: 3x interpolated in ln(k/X) vs log Z between Exo-FMS 1x and 10x, validated on Sonora 2020 at 1 dex
  (grey medians 1-3 %); hiT2 extension as for 10x.
- EOS composition eos_xh 0.7188, eos_yhe 0.2420, met/rad_met 0.4771; grav 1172.7924, x1 1.124210e10..1.595088e10.
- Solver/physics keys as 1x/10x; lhllc; ck_impl_conserve = 1.

## Steps
1. `git fetch fork`; build rt-integration **>= 2fd94098** (viper built 2fd94098; the later commits are docs only) for
   dhj (PROBLEM=deep_hot_jupiter_rt) with your incremental build_inc.sh (CUDA, HOPPER90). Record md5.
2. Package (2.8 MB): `git fetch fork data-w121-3x-pkg && git show fork/data-w121-3x-pkg:w121_3x_pkg.tgz > w121_3x_pkg.tgz`
   (md5 9dc9472ee4545b1b0150767577c068ec); `tar xzf w121_3x_pkg.tgz`; `bash w121_3x_pkg/SETUP.sh <repo root>`
   (links cia/ray/sw_flux from data/exo_fms_ck and writes w121prod_3x.athinput with absolute paths).
3. Smoke (50-100 cycles, as your w121prod_0928 smokes): rc 0, no FATAL/nan, ck 0 NOT-CONVERGED, dt ~15 s
   (viper: 15.19 s on 2 MI300A at nx1 76).
4. Production: 1 node x 2 H200 (as your 1x/10x w121prod), fresh start `-i w121prod_3x.athinput`, time/tlim = 300 rot
   = 3.304605e7 s, rst every 0.5 rot, same outputs; chain as you did for 1x/10x.
5. NOTE-2026-09-30-caltech-w121-3x.md: start estimate, job ids, binary md5, smoke; later rot-300 completion.
