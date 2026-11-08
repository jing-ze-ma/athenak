# NOTE 2026-11-08 Orion -> viper: HG MESA L_in(t) TASK not done (USER decision)

Re: TASK-2026-11-08-orion-hg-mesa-lin.md (f59a38312).

The USER decided that Orion will not do this task. Nothing was run and nothing is pending on Orion; the scratch
directory used for setup has been deleted.

Facts found during setup, in case the task is picked up elsewhere:
- Orion has no MESA install or MESA SDK (none under /orion/u/jinma, /u/jinma, /orion/ptmp/jinma; no module).
- The start model `2.754M_from_14.5M.mod` (or any S. Wu photo/.mod) is not on Orion either.
- MESA r15140 has `relax_initial_composition` and `relax_initial_entropy` (with `get_entropy_for_relax_from_eos`,
  `set_initial_age`) in star/defaults/star_job.defaults, and profile14 contains q, entropy, logT, logRho and all
  22 isotopes of approx21_cr60_plus_co56. Starting from a relaxed profile14 was considered and not pursued.

The consfix-1109 ledger TASK (4c95c5f43) is being run on Orion; its reply follows separately.
