---
name: cs-wb-arm-and-ke-gap
description: 2026-09-12 cs_mhd_prod3_wb (cs MHD + sp's WB block) DIED at rot 12.3 (dt collapse in 100 cycles, no precursor, no cell diagnostic on the MHD path); up to rot 12 it tracked plain cs_mhd_prod3 on KE/ME/1-ME/E to 10%, so WB is NOT what separates cs from sp. The "cs KE 7x below sp" gap is on the SP MHD side, not cs -- every cs arm AND sp HYDRO sit at KE 2-3e34 at rot 12, only sp MHD is 1.6e35 with rising total E
metadata:
  type: project
---

KE at rot 12 (dhj.mhd.hst cols 8-10; P = 3.05e5 s): cs_mhd_prod3 1.82e34, cs_mhd_prod3_wb 1.83e34,
cs_mhd_prod2 (WB + rot_potential) 2.13e34, cs_hyd_rs hllc/ausmpup 2.6-2.9e34, sp_dhj_hyd 2.34e34,
**sp_mhd_prod3 1.60e35** (plateau 1.2-1.6e35 to rot 150, total E climbing +5% by rot 8 then flat
near 7.35e38, while every cs arm and sp hydro DRIFT DOWN ~0.05%/rot).
So the anomaly is sp MHD = [[sp-mhd-energy-excess]] (cause open), not a cs deficit; the cs
comparison vs sp must use the sp HYDRO run or wait for the sp MHD excess to be understood.

WB arm death: bench/cs_mhd_prod3_wb, link 11618999, FATAL dt < 1e-2 at cycle 311578, t 3.741e6
(rot 12.27); dt 13-17 s until cycle 311300, 0.3 s at 311400. Floors normal in the last dhj.log row.
rst 24 (rot 12.0) and bin 6 (rot 12) exist for a restart reproduction. The mesh.cpp dt-collapse cell
diagnostic exists ONLY for hydro (ph->dt_diag); MHD has none -> add one before hunting the cell.
cs_mhd_prod3 (no WB) passed rot 23.6 clean at the same time.
