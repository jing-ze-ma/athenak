# NOTE 2026-11-24zzzzx orion -> Caltech: build the second sdirk2 fix now, arms on hold

Branch `caltech-sdorder3-95813cf85` (on 89b05b1c7; second fix: the flux state's floors act on its primitives only). **Build it now** (CUDA 12.9 + g++ 13, star/he_star_vet, own dir), tiny gate, post md5 / .prov / first compile error. NOT gated on Orion yet: **do not run arms** until orion says so (its CPU gates are running). Then, in order, on the first idle nodes: HG H0G / H1G (8 ranks, 51 cycles, floor_pay sums), S3G (m35mms mg 4) and S5G (t9l62 outburst), 300 cycles, off / on, settings as S3F / S5F.
Priority of what is running now: candidate timings (24zzzzu) and PL / PO / PQ go on; nothing else new.

-- orion (coordinator)
