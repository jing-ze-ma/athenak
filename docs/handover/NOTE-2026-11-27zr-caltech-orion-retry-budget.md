# NOTE 2026-11-27zr caltech -> orion: clean_retry_budget stored 0 FATALs every production restart on caltech-next / starmix

## clean_retry_budget: every production rst FATALs at restart on caltech-next-104c68b7e and caltech-starmix-bd1503bef

- **What happens.** A restart from any of our production rsts (bsg 00009, m35mms 00011, m35zams 00011, t9l62 00012) on the MIX binary stops at start-up on all 4 ranks:
  `### FATAL ERROR in src/rad_vet/rad_vet_args.hpp`
  `<rad_vet/clean_retry_budget> = 0: the key was removed, the code runs true (name clean_retry_budget=true to continue so)`
  This hit XR-bsg and all 8 MIX RST runs of 27zl.
- **Why.**
  - The production rsts store `rad_vet/clean_retry_budget = 0`.
  - In aa2646b92, rad_vet_params.cpp line 115 lists the removed key with its only value `=false`. With that rule, the stored 0 restarts fine on REL: RR-bsg and RR-m35 ran.
  - The auditfix retry-trap commit (afef77099 on caltech-next, 00052ca3f on starmix) changes that value to `=true`.
  - Reconfatal's rule then turns the stored 0 into a FATAL.
- **Scope.** Fresh starts are not affected. caltech-next-104c68b7e contains afef77099, so it has the same FATAL for every production restart; I expect it on the 27zq binary too (not run here, its suite has no production rst).
- **Workaround in my arms.** The MIX restart arms (XR-*, RST-MIX-*) name `rad_vet/clean_retry_budget=true`. That is the code's only behaviour, so the key moves no result. The REL arms (RR-*) name nothing.
- **For you to decide before the merge.** For example: accept a stored false/0 as the removed key's old meaning (no FATAL on restart), or keep the FATAL and add a release-note / RUN.md line that every production restart deck must name `clean_retry_budget = true`.
