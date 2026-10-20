# NOTE 2026-10-20 raven -> caltech: HOLD the fresh accretor (acc3d_trc_1020) — the deck changes to the production method (USER)

Audit (defaudit-1020): v11amwtr(+lz) pins an old method block copied from the 2-D accbig decks (time_scheme = be, implicit_gas_src =
matched, rad_inertia/nofb_work off, vet_bc_diff cell, vet_seg mean, f0_init zero, precond rbgs_fwd, ...). USER: drop it and use the
production defaults (second order in time, energy-conservative sources). Raven is building + screening deck v12lz now (hours).
- Do NOT start the fresh link 1 yet: hold / scontrol hold it if it is queued. Keep the OLD chain 5088058 running meanwhile.
- When v12lz is ready, a NOTE gives the deck (md5) — then start fresh on 3896d529 + v12lz. The race with Raven run12 stays.
- RSG E: unchanged (its fresh defaults are all on).
