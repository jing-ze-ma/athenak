# NOTE 2026-11-01 Caltech -> viper (USER): vet_gd twin bias check, HG and RSG E are CLEAN; only the accretor tip is affected

Follow-up to d07e7322.
- **Metric:** self-thick cells = |F_r|/(cE) < 0.01 (tau >~ 30), also checked at < 0.003. We compute D_rr - 1/3 and D_tt - 1/3 from the m1_vet dumps. Script csacc_1031/dthick.py; read-only.
- **HG** (prodG vet 48/43, t 240/215 d, ~2e7 thick cells):
  - p1-p99 within 1e-4 of 1/3;
  - |dev| > 0.02 in 0.08-0.09 % of cells, isolated, 1-12 % below r_ph under photospheric dips.
- **RSG E** (prodH vet 638/633, ~4.4e7 thick cells):
  - p1-p99 within 1e-4;
  - |dev| > 0.02 in 0.02 % (0.00 % at the strict cut).
- **Accretor sp wedge** (the sanity check that the metric catches the bug):
  - |dev| > 0.02 in 2.8 % of thick cells, systematically negative D_rr (p1 -0.043, min -0.075);
  - all in the L1-tip columns, 1.4-3.5 % below r_ph.
- **Conclusion:** HG/RSG need no restart fix. With geo42 + vet_rotate_avg 8 they are far less aspherical relative to the shell width than the tidal tip. The tip fix continues on csacc-1031 (twin0 / symref0 / twin_weight = j controls, then the thick-cell vanishing).
