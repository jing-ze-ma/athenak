# TASK 2026-11-17 viper -> Orion (USER): per-key CPU timing of the push-16 defaults

Thanks for af7917e12. Push 16 is 3-17 % slower than push 15 on CPU, and our rule is that a new default costing > 5 % stays opt-in.
Please bisect it on CPU.

- Head: fork clean-cc-push16 **5e33f056a**. It is 0aee9a208 plus your init-audit patch (4cacdae21, 8db6f6141) and your steady-start
  predict fix (5e33f056a).
- Same node and method as job 230771, decks flz_q, flz_Q, cs_closed, cs_feed1, cs_thin1, thin_f, RT feed, feed_x1, these arms:
  1. push 15 dfe21a746 (reference)
  2. 5e33f056a, defaults
  3. 5e33f056a with `clean_gas_exact = false`
  4. 5e33f056a with `clean_fv_joint = false`
  5. 5e33f056a with `clean_pc_factor = false` (matters only for x1line decks)
  6. 5e33f056a with all of them false (should match push 15 within noise)
- Also: is 5e33f056a with all new defaults off still bitwise to push 15 on the reg harness (A arm)? It should be, apart from decks that use
  steady start + predict.
- Reply as NOTE-2026-11-17-orion-viper-push16-perkey-timing.md.
