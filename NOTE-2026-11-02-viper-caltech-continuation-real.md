# NOTE 2026-11-02 viper -> Caltech (USER): the MHD continuation failure with roche_wb_tan is REAL (not a cmp artefact). Hydro is fine, as you showed.

Reply to f549c4948. Thanks for the hydro continuity test; the hydro launch is unaffected.
- **Re-checked with a header-independent comparison:** the bytes from `<par_end>` on (rstcmp2.py, Raven /raven/ptmp/jinma/merge_1101). Our earlier tail compare also used the LAST min(size)-100000 bytes, so it was not offset.
- The data sections have equal length, and ~half the bytes differ:
  - base d3f878be, r1m-type MHD deck (ohmic eos + STS), 150 + 150 vs 300:
    - DIFFER in every arm with roche_wb_tan = true (17.6M of 33.8M bytes);
    - IDENTICAL in the arms without it.
  - integ e32a104e and e1ecd210: 24.8M of 35.4M bytes differ.
- **So in MHD with resistivity, roche_wb_tan alone breaks restart continuation in the base code.** That fits your hydro PASS: the eta_b/ad_coef clobber that d0d1a738 fixes is MHD-resistive only.
- **Next:** our gate on 2a6b5e9fc (which has d0d1a738) and your MHD continuation pair will show whether that fix is all of it. Please use a `<par_end>`-based comparison.
- If it still differs, the next suspects are other MHD state the probe's ConToPrim touches and does not restore, e.g. the ck/radiation caches, wtemp warm starts, or C2P-dependent RKG register state.
