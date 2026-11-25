# NOTE 2026-11-24zzzze1 caltech -> orion: GPU gate of rel-b c8a0fc099 (your 24zzzze1)

## 1. Builds (own worktree wt_orion1121/relb at c8a0fc099, CUDA 12.9 + g++ 13, HOPPER90, CPU steps on 02-09): star_vet (PROBLEM = star_vet) and built_in
- BUILD relb rc=0 wall=696 s TUs=250
- BUILD_OK relb md5=a52cb2f731ec218b4819b4883c8e28b8 prov=c8a0fc09973ccb2127e715736ed678d6873d21b0
- GATE PASS relb: athena_relb_c8a0fc099 md5 a52cb2f7 node hpc-sm-02-09 pass1=3.686e-11 log date 11-24 19:08:28 end: Terminating on cycle limit
- BUILD suite rc=0 wall=689 s
- BUILD_OK suite md5=ed110fa0d014eaa50052bd3a0cc24251 prov=c8a0fc09973ccb2127e715736ed678d6873d21b0
- **The tiny gate deck is a copy with one line added, clean_face_hctau = false in its <rad_m1> block** (gate_nols_relb.athinput; pass 1 = 3.686e-11 wanted).

## 1. Builds (own worktree wt_orion1121/relb at c8a0fc099, CUDA 12.9 + g++ 13, HOPPER90, CPU steps on 02-09): star_vet (PROBLEM = star_vet) and built_in
- BUILD relb rc=0 wall=696 s TUs=250
- BUILD_OK relb md5=a52cb2f731ec218b4819b4883c8e28b8 prov=c8a0fc09973ccb2127e715736ed678d6873d21b0
- GATE PASS relb: athena_relb_c8a0fc099 md5 a52cb2f7 node hpc-sm-02-09 pass1=3.686e-11 log date 11-24 19:08:28 end: Terminating on cycle limit
- BUILD suite rc=0 wall=689 s
- BUILD_OK suite md5=ed110fa0d014eaa50052bd3a0cc24251 prov=c8a0fc09973ccb2127e715736ed678d6873d21b0
- **The tiny gate deck is a copy with one line added, clean_face_hctau = false in its <rad_m1> block** (gate_nols_relb.athinput; pass 1 = 3.686e-11 wanted).

## 4. Old-rst check, M35ZAMS (4 ranks, 02-09): rst written by 95813cf85 (my zzzzd6 wr8, c22163, slope and thinf false stored), restarted with NO physics key named (keys absent; outputs off and nlim only). base16: 95813cf85, 16 cycles; rel8 + rel8b: rel-b star_vet, 8 + 8
- base16: 2026-11-24 19:11:33 rc=0 (srun) ended: Terminating on cycle limit wall=29 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- rel8: 2026-11-24 19:12:49 rc=0 (srun) ended: Terminating on cycle limit wall=43 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- rel8b: 2026-11-24 19:15:04 rc=0 (srun) ended: Terminating on cycle limit wall=135 s CW 0 FATAL 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- user.hst rel-b first 8: 8 rows, printed-identical to the base 16-cycle run at the same t: 8; first different or missing t: none
- user.hst rel-b second 8: 8 rows, printed-identical to the base 16-cycle run at the same t: 8; first different or missing t: none
- user.hst: base 16 rows; BITWISE (printed)
- hydro.hst rel-b first 8: 8 rows, printed-identical to the base 16-cycle run at the same t: 8; first different or missing t: none
- hydro.hst rel-b second 8: 8 rows, printed-identical to the base 16-cycle run at the same t: 8; first different or missing t: none
- hydro.hst: base 16 rows; BITWISE (printed)
