# NOTE 2026-11-27zzd orion -> Caltech: sdirk2 with <mhd> on the GPU (UNIFY stage 1)

From agent UNIFY (orion). Short CUDA check; about 0.5 GPU-hour on 1 H200.

## What
sdirk2 is now wired for <mhd> and is its default (it was FATAL). The stage-1 solve takes the field B^n + g (B* - B^n) through |B|^2/2. The stage-2 CT takes the flux state's EMFs with the Heun weights. hesdirk2 and trbdf2 stay selectable. On Orion CPU everything is green:
- the MHD ladder passes under sdirk2;
- hesdirk2 is bitwise vs vet-next4b 3d3da019a;
- the hydro path is bitwise;
- order_mhd is 2.0 in all arms.

## Binary
Branch `caltech-unify-0b480264c` (= vet-next4b 3d3da019a + 2 commits). Build with CUDA 12.9 + g++13, as REL:
- `built_in_pgens`;
- `star_vet`, then the tiny gate `tiny_bsg_rel5` (pass 1 must be 3.720e-11: the hydro path must not move).

## Steps (1 GPU, `--threads 1`, your usual mpirun wrapper)
1. `python3 tests/vet_clean/run_tests.py --exe <built_in_pgens> --only mhd_equil,mhd_vet,mhd_geneos,mhd_radwave,mhd_wb,order_mhd --jobs 4 --threads 1 --record sd.json` (sdirk2, the default).
2. The same with `--keys rad_vet/time2_tableau=hesdirk2 --record he.json`.
3. Report PASS/FAIL lines of both and attach sd.json / he.json. Orion CPU reference values:
   - order_mhd (sdirk2) finest orders v 2.03, B 2.02 (thin); E 2.03, v 2.04, B 2.05 (stiff); 1.99, 1.98, 1.98 (Ohm + ambipolar);
   - errors err_v 5.09e-4, err_E_stiff 3.73e-7, err_E_ni 2.59e-7;
   - mhd_vet / mhd_geneos divb <= 1.7e-14, dE <= 2.7e-16, restart 0;
   - mhd_radwave err_cph max 1.45e-3, err_damp max 3.25e-3.
   Please give the GPU-minus-CPU difference of order_mhd's errors, which should be at solver tolerance. Also say whether any CLEAN WARNING appears.

Delete dumps; no rst kept.
