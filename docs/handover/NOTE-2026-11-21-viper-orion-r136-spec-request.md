# NOTE 2026-11-21 viper -> orion (USER request): please read the R136a1 ic_star spec on Caltech — which `mesa_profile.top`?

Follow-up to the T9 ic_star thread (NOTE ...-ic-builder-t9-answer, 8d66a59df). The USER asks which top mode the R136a1 production
column used. The column is ic_r136a1_B.txt (md5 ec3c8a81), source of viper's E2 IC ic_r136a1_B_dfl17.txt (e8de8c76). It was built
on Caltech in `/resnick/groups/carnegie_poc/jingze/r136a1_1104/`, and its spec is not in the repo. The IC header says only
"tools/ic_star r136a1 column, seam fixed point (r136a1_1104/tools/seam_only.py) ... seam pass 2: rows = the code dump
.../r136a1_1104/ic/BB/dump_2.txt (sph1d re_atm hse), scaffold of the dump (taper off)".

Please report, read-only:
1. The spec file(s) used for that column (path + full text), especially `mesa_profile.top` (hopf / surface; default hopf),
   `mesa_profile.top_fit`, `mesa_profile.r_max_rsun`, `taper`, `seam`, `startup.re_zone`, `grid.tauj`, opacity block.
2. The make_star_ic.py log of that build: the line "mode mesa: tau-3 join at r_j ..." (or "(surface)"), the Hopf-layer
   result (R_ph, rho_ph), and any FLAG lines (Gamma, branch).
3. seam_only.py (path + text or a summary of what it changed vs the tool's own seam loop), and how pass 1 -> pass 2 went.
4. If the R136 column did go through the outward Hopf layer to tau 2/3 despite its Gamma > 1 layer at 25.0-25.3 Rsun, the
   Gamma along the Hopf layer (or the r_j / R_ph values), so we can compare with T9's failure at r_j = 83 Rsun.

Purpose: decide whether T9L6.2 should use `top: surface` + r_max cut, or whether R136 shows another route. No runs needed.

-- viper (for the USER)
