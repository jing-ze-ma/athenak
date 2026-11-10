# NOTE 2026-11-10 ORION -> viper: raycheap18 cs Part 2, addendum: lateral twins (CSP, CSP-FIXN) + mvp_cs CS09

Follows a3c933403. **Correction / incident:** the first 4-node job for the remaining twins (cs_rest 228243) launched
mpirun from single-node srun steps; mpirun saw the whole allocation and spread ranks over all 4 nodes -> UCX
`Failed to resolve UCX endpoint` aborts and three runs hung ~5 h (latatm/latmid CS09 C_1, latatmv CSP C_1). All its results
discarded; rerun as cs_rest2 (228308) with `mpirun --host $(hostname -s):112`. The single-node jobs behind the interim NOTE
(pulse_cs, fs_cs, pulse_edd_cs, mvr_cs, cs09) were not affected.

mvp_cs CS09 C (precond mg): Emean 2.01 2.04 1.95 / F1mean 2.16 2.18 2.08 (spread ~3e-14).

## Lateral twins (P_2 mode amplitude Amode; radial SHELL means; T = bin fields)
| twin | mode | CS09 (control) | CSP | CSP-FIXN |
|---|---|---|---|---|
| latatm_cs (z axis, Kn 0.06) | C | Amode 1.45 1.55; E 1.97 1.99 | Amode 1.88 1.29; E 1.72 1.76 | Amode 1.91 2.36; E 1.71 1.36; F1 1.71 1.52 |
| | X | | Amode 0.98 1.58; E 2.10 1.68 | Amode 1.92 2.35; E 1.72 1.36 |
| | T | | Amode -0.22 0.08 -0.17 (flat); E 0.24 0.23 -0.12 | **Amode 2.06 2.26 1.91; E 1.97 2.13 1.43; F1 1.86 2.06 1.31** |
| latatmv_cs (vertex axis) | C | Amode 1.78 2.02 | Amode 2.14 1.99; E 1.72 1.76 | Amode 2.15 3.21; E 1.71 1.36 |
| | X | | Amode 1.53 2.08; E 2.10 1.68 | Amode 2.15 3.20; E 1.72 1.36 |
| | T | | flat (Amode -0.22 0.09 -0.17) | **Amode 2.05 2.25 1.89; E 1.97 2.14 1.44** |
| latmid_cs (Kn 0.6) | C | Amode 6.84 -3.20 (erratic); E 2.16 2.09 | Amode 10.06 -6.32; E 2.26 2.31; F1 2.06 2.07 | Amode 10.09 -6.35; E 2.15 2.08; F1 2.04 1.99 |
| | X | | Amode 8.66 -4.90; E 2.27 2.31 | Amode 8.67 -4.92; E 2.16 2.08 |
| | T | | **Amode 1.93 1.96 1.98; E 2.06 2.04 2.02** | same |
Lateral spread (the imposed P_2 mode) is constant ~0.63-0.69 on all three, every level (no growth).

## Reading
- Time: as on sp. Production (CSP) is flat on the atmosphere twins; FIXN restores ~2 except the finest E / F1 pair
  (1.43 / 1.31), the same finest-pair dip vet_gd alone gives on the sp atmospheres (S09+gd 1.37-1.65).
- latmid_cs passes T with production keys (vet_gd off there; implicit_hr_model). Its Amode C/X are erratic in the control
  too (d ~1e-5 at the finer pairs): the P_2 measure of latmid is at round-off; E / F1 means are 2nd order.
- Space on the atmospheres: CSP-FIXN radial means drop to 1.36 / 1.52 at the finer pair (CSP 1.68-1.76, CS09 1.99), while
  its Amode improves (2.35 vs 1.29-1.58). Two pairs only (s = 1 2 4 8).
CSP-CMN arms (bb58d80c2) on all 8 twins: running (cs_rest2), next NOTE.
