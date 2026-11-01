# NOTE 2026-11-02 viper -> Caltech (USER): with your face-EMF fix, the MHD restart continuation PASSES on our side (integ-test-1101 = af982e409)

- **integ-test-1101 = af982e409** = 2a6b5e9fc + cherry-pick of your 267a32d0 (the probe restores the six face-EMF arrays).
  - The only conflict was with our xe_init_pending restore line; we kept both.
- **Raven CPU gates on af982e409** (r1m-type MHD deck, ohmic eos + STS, roche_wb_tan and all defaults on):
  - (d) 150 + 150 vs 300 cycles: rst data IDENTICAL, both in the tail compare and in the `<par_end>` payload compare;
  - (a) old restart without keys vs d3f878be: IDENTICAL;
  - (b) fresh no-key vs all-explicit: IDENTICAL;
  - (c) explicit old values vs base: IDENTICAL.
- **One open detail:** the hst files of the continued run still differ from the straight run in the second half, although the final state is bitwise identical. We are checking whether this is hst bookkeeping (row timing or a duplicated row at the restart) or a real diagnostic difference. We will NOTE.
- **Please send:** your GPU continuation pair on 267a32d0 and the probe-audit result. If both are clean, we push af982e409 to fin51-r-fix. That is your r1m launch build.
- **Note: the base fin51-r-fix d3f878be still has both probe bugs** (eta and face EMFs), so any MHD run with roche_wb_tan + resistivity on d3f878be is non-continuous at restarts and spikes on a fresh start with the x_e tracer + ambipolar. Hydro is unaffected.
