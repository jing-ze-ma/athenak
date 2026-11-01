# NOTE 2026-11-01 viper -> Caltech (USER): final integration commit e1ecd2101 (fork/integ-test-1101)

- **integ-test-1101 = e1ecd2101** = e32a104ea + merge of fin51-r-fix d3f878be7 (31d5c611a) + the relax_drag = false fresh-start default.
- **One merge conflict, in the dhj flux-hst labels.** Column order is now ... Ttide (roche), Tdrag (am_hst), then Mdrain/Edrain (dropped with a NOTE if > 24). With am_hst on, Mdrain/Edrain sit one column further right than in e32a104e. State is unaffected.
- **Please run on e1ecd2101 (short):**
  - (a2) key-free bitwise on the same developed r1m restart vs your d3f878be binary, ~150 cycles: rst data + dhj.mhd.hst identical;
  - optionally a fresh no-key r1m start, ~50 cycles, to confirm the defaults block now also shows relax_drag = false.
- Your (a)/(b) on e32a104e stand. On PASS we push e1ecd2101 (or its successor if CPU gates force a fix) to fin51-r-fix.
