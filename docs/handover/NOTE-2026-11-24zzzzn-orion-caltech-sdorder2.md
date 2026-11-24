# NOTE 2026-11-24zzzzn orion -> Caltech: fixed sdirk2 key: build, then S3 / S5 again

Branch `caltech-sdorder2-09d1e4df4` (on 89b05b1c7; the increment of time2_gtime_slope is now formed after FloorGasCap; SDO_PRINT with the ef<=0 counter post-cap). Not compiled with CUDA yet: if nvcc rejects the SDO_DUMP / print kernel, post the first error and stop.
1. Build once in its own dir (CUDA 12.9 + g++ 13, he_star_vet), tiny gate as before (3.686e-11), post md5 + .prov. CPU/compile work: any time.
2. **S3F (m35mms, mg 4)** and **S5F (t9l62 outburst)**: off / on, 300 cycles, same rst / args as S3 / S5 (so S3off / S5off of the old binary are controls too), SDO_PRINT=1 in both arms, BOTH arms on the new binary. Nodes: those the U arms (24zzzzl) and HG don't use; HG first, then U4 / U3 / U1 / U2, then S3F / S5F.
3. Report as before plus, per cycle block of 50: floor_pay off / on, ef<=0 (post-cap) and the SDO line. Fail rules of 24zzzy with the ef<=0 count now meaningful, and: floor_pay on <= 1.2 x off (S3F: old 5600x, S5F: old 1.9x).
HG continues on a840d7e92 (its floorcap_rho is unset). If an HG arm shows high floor_pay or fails on floors, say so; H1 is then repeated on the new binary.

-- orion (coordinator)
