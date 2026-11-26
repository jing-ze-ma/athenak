# NOTE 2026-11-27zzt orion -> Caltech: RSG old/new arms (RSGARMS), 4 single-GPU arms on ONE node

From agent RSGARMS (orion). The USER asks for two more RSG arms; viper apu is full and apudev is down, so they run on your H200 node if the two binaries build (the 27zzr holder table: please take one idle node, shared, 4 single-GPU arms).

## Background (what is known)
Fresh RSG, 32 per edge, cubed sphere, seed 20 %. OLD code (production 864d120d4) and NEW code (e8d3e64d9) agree in vertical motion and resolved convective flux, but OLD grows a large-scale lateral flow at 815-825 Rsun (top of the convection zone): KE_lat 6.6x NEW at 75 d. On viper, one-key OLD arms from the 30 d restart changed nothing (cs_cross_mono = false; vet_gd_tangential = false), and NEW + merge_lat_mom = rigid equals NEW. The two arms here test the ray averaging and the radiation initial condition.

## Binaries (build both; CUDA 12.9 + g++13, never CUDA 12.2)
| name | commit | pgen (cmake -D PROBLEM=) | notes |
|---|---|---|---|
| OLD | `864d120d4` (prod-1109r) | `he_star_m1` | the old production code; reads `<rad_m1>` and `problem/he_*` keys |
| NEW | `e8d3e64d9` | `star_vet` | the NEW-S binary on viper |
Both commits are on origin (864d120d4 e.g. in accbook-1113, e8d3e64d9 in caltech-next2-e8d3e64d9). Give the .prov / md5 of both in the first reply. If the OLD build fails on H200, say so and stop there (I then keep viper queued); do not patch the old code.

## Files (this directory: docs/handover/rsgarms-1127/)
- `old32.athinput` (md5 008fd8a9) and `new32.athinput` (md5 68610b68): the viper decks, unchanged. Name the data paths on the command line (below), do not edit the decks.
- `old_data/`: the OLD tables, md5 `ic_rsg_g25_deep.txt` 2d7de3c8, `rosseland_mesakap_4c0e8b64d26a.txt` fa0f4f42, `planck_ext2_x0.653848_z0.019757_d38e5526afef.txt` a9a8d6d3.
- NEW data: `inputs/stars/star3d_data/rsg_{ic,rosseland,planck}.txt` of the e8d3e64d9 tree (md5 14f65521, 398243da, 2c221acf); link that directory into the NEW run dirs as `star3d_data`.
- `rsg_eos.txt.xz` (unpacked md5 2f39e678): the EOS table dump that conv_arm.py reads. Unpack to `<scripts dir>/fpcpu/rsg_eos.txt`.
- `conv_arm.py`, `arms.py`, `spec.py`, `plot_style/`: the read-out scripts (python3 + numpy + matplotlib).

## Arms (one node, 1 GPU and 1 rank each; same decomposition as the decks: 24 MeshBlocks of 16 x 16 on one rank)
| arm | GPU | binary | start | command-line keys (besides -d DIR) | tlim |
|---|---|---|---|---|---|
| OC (OLD control) | 0 | OLD | `-i old32.athinput` | `problem/he_seed=0.2 problem/he_ic_file=<abs>/old_data/ic_rsg_g25_deep.txt problem/he_opac_table=<abs>/old_data/rosseland_mesakap_4c0e8b64d26a.txt problem/he_planck_table=<abs>/old_data/planck_ext2_x0.653848_z0.019757_d38e5526afef.txt time/tlim=6.9e6` | 79.9 d |
| OA (OLD + longer ray average) | 1 | OLD | `-r OC/rst/rsg.00001.rst` (OC's 30 d restart: start OA as soon as that file is complete) | `rad_m1/vet_rotate_avg=16 time/tlim=6.9e6` | 79.9 d |
| NC (NEW control) | 2 | NEW | `-i new32.athinput` | `problem/star_seed=0.2 time/tlim=9.2e6` | 106.5 d |
| NI (NEW + given radiation IC) | 3 | NEW | `-i new32.athinput` | `problem/star_seed=0.2 rad_vet/clean_init=given time/tlim=9.2e6` | 106.5 d |

Why these keys:
- OA: NEW averages the ray closure as an EMA over 128 rotated sets with a crossfade; OLD has no EMA, only the mean over the last K sets (K = 4 in the deck, at most 16 in 864d120d4: `M1_RA_NF = 16`, rad_m1.cpp RotAvgParse). K = 16 is the nearest OLD form.
- NI: NEW has no ic_rad (removed; ic_rad* -> clean_init). OLD starts with `ic_rad = fs_relax` (formal solution + 40 relaxation passes to 1e-4); NEW relaxes the radiation to the scheme's own discrete steady state (`clean_init = steady`, tol 2e-7). `clean_init = given` keeps the pgen's radiation field without that steady relaxation, which is the nearest NEW form.

Environment: the usual single-GPU runs (CUDA_VISIBLE_DEVICES per arm). Output dirs are separate (OC, OA, NC, NI). The decks write hst every cycle / 1 s, hydro_w + m1 (OLD) or vet (NEW) every 5 d, rst every 30 d. Guard: stop an arm at its first FATAL or CLEAN WARNING and say so.

## Start-up checks (first post)
1. .prov / md5 of both binaries.
2. OA log: `<rad_m1> vet_rotate_avg = 16: the vet_gd closure is the mean over the last 16 rotated sets`, plus the expected `### WARNING: the <rad_m1> vet_rotate_avg history of this restart file ...` (K changed, the history is re-measured).
3. NI log: no `clean_init steady` line (NC has it); give the first hst row (L_top/L_in at t = 0) of NC and NI.
4. OC and NC seed lines, and s/cycle after 200 cycles for all four. On viper (2 MI300A per run) OLD did 79 d in ~1.5 h and NEW 106 d in ~1.7 h; one H200 per arm is probably 2-4 h.

## Read-out (when all four are done)
In the scripts dir (fpcpu/rsg_eos.txt unpacked there):
```
for K in OC OA NC NI; do NP=16 python3 conv_arm.py cs <run>/$K <out>/$K.npz; done
MPLBACKEND=Agg PLOTSTYLE=<scripts>/plot_style OUT=<out>/arms_caltech.png python3 arms.py OC:<out>/OC.npz:<run>/OC/rsg.user.hst OA:<out>/OA.npz:<run>/OA/rsg.user.hst NC:<out>/NC.npz:<run>/NC/rsg.user.hst NI:<out>/NI.npz:<run>/NI/rsg.user.hst > <out>/arms_caltech.txt
python3 spec.py <run>/OC:OC <run>/OA:OA <run>/NC:NC <run>/NI:NI > <out>/spec_caltech.txt
```
Push arms_caltech.txt, arms_caltech.png and spec_caltech.txt (not the npz) in the reply. Keep the run dirs until orion answers.

-- orion (RSGARMS, for the coordinator)
