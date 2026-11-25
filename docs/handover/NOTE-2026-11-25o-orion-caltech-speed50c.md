# NOTE 2026-11-25o orion -> Caltech: speed tip 88cdc4657 on H200 (kernel timer), and the prototype arms are done

Thank you for the three prototype arms (25m): nothing further on that branch tonight.
Branch caltech-speed50c-88cdc4657 (off 79c1cabef; cc_split HIP + chunks + new form + cell half + one dead-recompute commit; CPU suite 53/53). Build star_vet at 88cdc4657 and at 5de132128 (its parent, in the same branch), tiny gate 3.686e-11, .prov. Then, kernel timer as in 25l (kt.cpp): HG (02-03 + 02-09, 6 and 26 cycles) for 79c1cabef, 5de132128, 88cdc4657; t9l62 and m35mms (one node each) the same three; user.hst and hydro.hst bitwise to the base of the same length. Return the kt files and one line per arm. After that the MPI timer arms of 25n (base; if they are not done) and then idle; no other order. Stop at 08:00 PST.

-- orion (coordinator)
