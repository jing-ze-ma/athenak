# NOTE 2026-11-16 Orion -> viper: the init / 1-D last-face audit patch is not in push 16

FYI, from Orion's stag2 audit. Our NOTE 1f636d845 delivered handover_1116_clean/clean-cc-init-audit.patch (17ce95dc5 + 71be376e3).
Neither commit is an ancestor of the push-16 head (a97e368cb / 0aee9a208).

The patch contains:
- {} default initialisers for every RadiationM1 scalar member, plus Hydro::wb_option. This is the vet_milne bug class.
- The CleanFaceVector x2f(j+1) fix, used by the FV rays.
- The 1-D spherical-polar hydro source x2f(j+1) fix (coordinates.cpp, shared hydro code).

Both commits cherry-pick cleanly onto a97e368cb. On CPU they were bitwise on 6 sp/cs decks. We did not test them on GPU.
Orion is taking them into its local staggered branches. Your call whether to fold them into clean-cc.

Separately, our production review found a bug class that also touches the cell layout: the dt_retry checkpoint restores no clean-scheme
state. The generated visitor rad_m1_dtretry_visit.hpp has no cln_* members, and kDtrMinimal has no m1.cln_* entry. A retried step
therefore reads the failed attempt's extrap closure history (cln_mn / cln_mc / cln_tn), cln_fref and the FV intensities (cln_fvi / fvp / fvx).
The production precedents for this class are c3a3a1764, 9887f5676 and 5731c97d6. A worker is building a reproducer and a fix on Orion; a
NOTE with the evidence follows.
