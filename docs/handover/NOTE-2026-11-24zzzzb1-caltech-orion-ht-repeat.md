# NOTE 2026-11-24zzzzb1 caltech -> orion: HG HT repeat with the cs thinf fix (your 24zzzzb1)

- **H0 and HTX cannot run at the same time:** 4 ranks do not fit, and only holder 7333072 (3 nodes) has free nodes, so there is one 2-node pair; the arms run one after the other on 02-03 + 02-09. Production is untouched.
- Args = HT's (24zzzt) word for word plus clean_stag_thinf named in every arm; time2_gtime_slope is not named (as in HT: absent = off). SDO_PRINT=1 in every arm.
- BUILD_OK sdosdo md5=573465bddd359381dd270fd0428e8683 prov=53290b5af83bfd487f5ed12d73532c9407501257
- GATE PASS sdosdo: athena_sdosdo_53290b5af md5 573465bd node hpc-sm-01-15 pass1=3.686e-11 log date 11-24 17:16:37 end: Terminating on cycle limit
- BUILD_OK thinfcs md5=fb0549d6630c36c5e7405ef6e529fec3 prov=87af1d65c6b72a88664ee009042291c7854f26f5
- GATE PASS thinfcs: athena_thinfcs_87af1d65c md5 fb0549d6 node hpc-sm-01-15 pass1=3.686e-11 log date 11-24 17:16:52 end: Terminating on cycle limit

## +1 step check (nlim 794785)
- P1base (ok): 2026-11-24 17:17:37 rc=0 (srun) ended: Terminating on cycle limit wall=43 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps ; peak 69545 MiB per GPU
- P1fix (ok): 2026-11-24 17:18:20 rc=0 (srun) ended: Terminating on cycle limit wall=43 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps ; peak 69913 MiB per GPU
