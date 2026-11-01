# w121_s10fresh_1027 package: fresh-start magnetized WASP-121b (twin of Raven arm w121_mhd_1003/s10_e15)

Prepared on viper 2026-10-27. Verify the files with `md5sum -c MD5SUMS` before you start.

## Contents
| path | what | read by |
|---|---|---|
| `s10fresh.athinput` | the deck. Fresh start t = 0, C32 x nx1 112, x1max 2.1e10, 3x metals, 10 G dipole, ohmic eos resistivity (max_eta 1e15), RKG2 STS, lhlld | `-i` |
| `ic/ic_w121_3x_f5.txt` | 3x face-5 RCE IC, T(p) (the same file w121fresh_3x_f6 uses) | `problem/ic_profile` |
| `ck/CK3X/` | ck tables: `ck/Premixed_3x_g8_11_hiT2.txt`, `CE_tables/FastChem_ck_3x_int_hiT2.txt`, `sw_flux/sw_band_flux_W121_11.txt`, `cia/`, `ray/` | `problem/ck_table`, `ck_data_dir`, `ck_ce_table`, `ck_swflux` |

The general EOS (`general_eos = table`) is computed at start-up from the deck's eos_* keys. It reads no file.

The table paths in the deck are **relative**: `ck/CK3X/...` and `ic/...`. Start athena from a run directory
that contains `ck` and `ic`, either as copies of these directories or as symlinks to them. Do not use `-d`.

## Build (Caltech H200)
- Source: `fork/fin51-r-fix` = **6fac5c7ff** (the binary Caltech already runs for he_star and the accretor).
- Build it incrementally in your existing nofma tree, as for those two runs:
  - `PROBLEM=deep_hot_jupiter_rt`
  - `Kokkos_ENABLE_CUDA=On`, `Kokkos_ARCH_HOPPER90=On`, `Athena_ENABLE_MPI=ON`, Release
  - CUDA 12.9, gcc 13.2, hpcx, host `-ffp-contract=off` (nofma). Use the same flags as the he_star_m1 6fac5c7f build; only PROBLEM differs.
  - Name: e.g. `athena_gpu_deep_hot_jupiter_rt_6fac5c7f_nofma`.
- viper reference binary (MI300A, ROCm 7.2): `builds/bin/athena_dhj_gpu72_6fac5c7f` md5 8bafe0d918cccfc4cd202555bc8ddae5.

## Run (1 node, 4 H200, 4 ranks; 24 MeshBlocks = 6 per GPU)
```bash
cd <rundir> || exit 1          # contains ck -> <pkg>/ck and ic -> <pkg>/ic
# smoke: about 200 cycles
srun -n 4 <athena> -i <pkg>/s10fresh.athinput time/nlim=200 > smoke.log 2>&1
# production (fresh run dir, no overrides): stops at the wall limit with a final rst
srun -n 4 <athena> -i <pkg>/s10fresh.athinput -t <hh:mm:ss = holder end - 20 min> > run.log 2>&1
# continue: srun -n 4 <athena> -r dhj.NNNNN.rst -t ...
```
Outputs:
- `dhj.mhd.hst`: every 0.01 rot = 1101.54 s.
- `dhj.user.hst`: the flux, wall, floor and angular-momentum columns.
- `dhj.log`: every 0.1 rot.
- bins `mhd_w_bcc`: every 0.5 rot.
- rst: every 1 rot = 110154 s.
