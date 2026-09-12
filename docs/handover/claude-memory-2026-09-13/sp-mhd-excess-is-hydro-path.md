---
name: sp-mhd-excess-is-hydro-path
description: 2026-09-12 bench/sp_excess (6 arms x 0.5 rot, sp_mhd_prod3 input, binary 9a9396f7): the sp MHD energy/KE excess NEEDS NO FIELD -- mhd_b0 (bbot=0, ME=0, tfloor 177) reproduces mhd_ctl to 1% (E-E_hyd +4.0e36 at 0.5 rot, horizontal KE 7.6x hydro) => it is the MHD MODULE'S HYDRO PATH vs the hydro module; uniform theta cuts it 6.5x (still 1.46x KE); low-beta fallback bitwise inert; mhd_b0 dt 16 s vs hyd_ctl 6.9 s
metadata:
  type: project
---

Arms (bench/sp_excess/<arm>/out, analysis/excess.py): mhd_ctl, hyd_ctl (same input, <mhd>-><hydro>,
hlld->hllc, MHD-only keys dropped), mhd_b0 (bbot=0), mhd_thuni/hyd_thuni (f_stretch_theta=1e-3),
mhd_lowbeta (cs_lowbeta_fallback=0.5: BITWISE = mhd_ctl, never engages or inert in 9a9396f7).
E - E_hyd at 0.1/0.25/0.5 rot: mhd_ctl +4.5e35/+1.5e36/+3.9e36, mhd_b0 +4.5e35/+1.5e36/+4.0e36,
mhd_thuni(vs hyd_thuni) +6.7e34/+2.5e35/+6.1e35. KE_h ratio: 2.1/4.8/7.6 (ctl and b0), 1.09/1.48/1.46 (thuni).
Mass: MHD arms +6e22 (gain), hydro arms -4.7e23 (loss), same IC 3.60049e26. dt at 0.5 rot: mhd_ctl 2.4,
mhd_b0 16.0, hyd_ctl 6.9, mhd_thuni 1.95, hyd_thuni 5.25.
Floors cumulative: mhd_b0 dfloor 6.6e8, efloor 1.6e7, tfloor 177 ~= hyd_ctl (7.2e8/1.7e7/4302): NOT floors.
**So:** every earlier "sp MHD excess" reading (temperature floors, nightside field 7x) was a symptom;
the cause is a difference between hydro/ and mhd/ module code paths on sp with B=0: candidates = the
Riemann solver (hlld at B=0 vs hllc), polar treatment (polar_emf_diss, polar averaging), WB source
(wellbalance in mhd_update vs hydro_update), c2p floor path (cs-mhd-c2p-floor-corrupts-ue), the pgen
source hooks branching on pmhd vs phydro, conduction/RT coupling through <mhd> vs <hydro>. The old
"(sp) MHD vs hydro" comparisons [[sp-hydro-vs-mhd-comparison]] [[sp-mhd-energy-excess]] measured this
module difference, not the field. cs MHD == cs hydro, so the difference is sp-specific (polar/stretch).
NEXT: dumps of mhd_b0 vs hyd_ctl at 0.25/0.5 rot (where is the excess: r, lat, day/night); short 0.1-rot
GPU arms flipping solver / polar / WB / conduction on BOTH modules.
