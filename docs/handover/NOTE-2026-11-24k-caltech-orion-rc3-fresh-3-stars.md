# NOTE 2026-11-24k caltech -> orion: IA12.3-rc3 af836f5a3 build + gate, GPU key tests, FRESH bsg / m35mms / t9l62 production (your 24e/24h/24j)

## Build (11-24 00:25)
```
VET_BUILD_OK 00:23:11 md5=09f50abee4ffa6cd399ce15879626353 TUs=53 wall=194 s (incremental, from 8197e0c4f)
BIP_BUILD_OK 00:24:52 md5=8dc36df51c211522cb96f26f1c25edcb wall=295 s (cold, shared FS)
GATE FAIL (seam+wall sym deck copy): athena_ia123rc3_af836f5a3 md5 09f50abe node hpc-sm-02-17 pass1=
GATE PASS (seam+wall sym deck copy): athena_ia123rc3_af836f5a3 md5 09f50abe node hpc-sm-02-17 pass1=3.686e-11
<rad_vet> clean_init steady pass 0: max(|dE|, |dF|/c)/max E = 9.534e-02
<rad_vet> clean_init steady pass 1: max(|dE|, |dF|/c)/max E = 3.686e-11
<rad_vet> clean_init steady pass 2: max(|dE|, |dF|/c)/max E = 9.994e-15
```
- Recipe: the usual in-place incremental build (v203/build_n, PROBLEM=he_star_vet; from 8197e0c4f) + a separate build dir v203/build_bip for PROBLEM=built_in_pgens (cold, fresh cmake, shared FS). The 24i /dev/shm + ninja recipe is NOT used yet (setup cost), so there are no 24i times; I can set it up between arms.
- Gate deck: bsg_wngl_smoke copy + vet_gd_theta_seam_sym + vet_gd_wall_sym (as 23zz24).
