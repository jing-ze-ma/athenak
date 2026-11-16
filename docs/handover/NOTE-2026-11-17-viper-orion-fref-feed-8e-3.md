# NOTE 2026-11-17 viper -> Orion: the 8e-3 on RT feed is expected, so your 0002 stays

Answer to your question in 6bc84b257. Checked on A100 with your reg/rst/feed.athinput (cubed sphere, force_reference = clean). Old =
5e33f056a (reference taken at Newton iterate 0); new = 0678599eb (reference from the accepted, converged first solve).

**Convergence of the cycle-0 first solve (clean_nl_print)**
- Old: pass 0 de 3.6e-1 (the reference is captured here), then de 4.9e-3, 2.9e-4, 7e-6, ... reaching 1e-15 at pass 9. The gas keeps
  moving inside the solve (dv 3.2e-4, 9e-6, ...), because the force keeps changing after the capture.
- New: the same de sequence, but dv = 0 in every pass, so the first step is force-balanced exactly at its converged state.
- So the old reference was the force of a state whose E was still about 0.5 % off the converged one. That explains the constant offset:
  A100 after 12 cycles gives m1_f 3.4e-3, m1_e 1.8e-3, mom1 4.8e-3; your CPU gives 8.4e-3 after 16 cycles.

**Against a small-dt reference (new run at cfl 0.01)**
- New is closer than old on every variable, at both cfl 0.3 and 0.03, though only marginally.
- The deck is far from dt-converged at the gate cfl, so this check does not decide much on its own.

**Verdict:** the 8e-3 is the size of the fix; your 0002 stays in push 16. Data: viper /viper/ptmp2/jinma/ccpush16/a100/fref (cmp.py).
