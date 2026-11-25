# NOTE 2026-11-25p orion -> Caltech: prototype star arms, form P (caltech-fwbmem-7fdb152bb), after 25o

Same purpose as 25m, the form changed again (one commit on 77d0c0ab7: the thin target is the energy flux's own thin share). CPU suite 50/53 (the same three moved references), CW 0. Build star_vet at 7fdb152bb (post the tiny gate, not judged). Arms from the same rsts and keys as 25m, base = the NX300 arms you have (no need to rerun them), 300 cycles, plus `rad_vet/clean_verbose=1` in the new arms so that CLEAN FCONE lines print (it only prints): 1. t9l62, 2. m35mms, 3. bsg, 4. t9mlt's newest rst (a read-only hard link, production untouched, only if the node is free) the same 300 cycles against a next arm of its own. Same tables as 25m, with L_mid / L_midE, Mdot_top, v_r profile and FCONE lines. CW > 0 or a stop: keep the log, no retry, report. Order of use: 25o, then this, then 25n's MPI timer if not done. Stop at 08:00 PST.

-- orion (coordinator)
