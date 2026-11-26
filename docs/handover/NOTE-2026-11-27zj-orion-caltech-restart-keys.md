# NOTE 2026-11-27zj orion -> Caltech: restart-key warning, check 4 stars on their production rsts (RECONFATAL, idle nodes, CPU is fine)

From agent RECONFATAL (orion). **Idle nodes only; never stop or share a node with production.** Each check is 3 restarts of a few cycles each.

## Why
USER decision (11-26): when a restart takes a key's value from the rst because the `-i` deck does not name it, the code now prints one line per key
`CLEAN WARNING restart: block/key = <stored> from the restart file (code default <d>; the deck does not name it)`
and only when the stored value differs from the code default. It is a counted CLEAN WARNING, never FATAL, and results do not change. A restart without `-i` prints none.
To keep production restarts at zero warnings, every deck must name those keys at the values the rst stores. On Orion the bsg, star3d_sp and star3d_cs templates printed 21 / 21 / 18 lines. With the keys named they print 0, and fresh runs, rst bodies and restarts are bitwise unchanged.

## Binary
Branch **`caltech-reconfatal-c9805a084`** (pushed with this NOTE) is aa2646b92 plus 2 commits:
- 4f9931fbd adds FATALs for keys that would not act (reconstruct other than plm on spherical polar or with wb_dirs = all; polar x2 bcs without use_polar_boundary; wb_x1 / wb_rho / wall_wb_onesided without wellbalance_dynamic; mhd fofc_rsolver = hlle without fofc; mhd polar_emf_diss without a polar boundary; dt_retry_views = all on the cubed sphere) and the restart warning.
- c9805a084 names those keys in inputs/stars.

Orion gates: tiny_bsg_rel5 PASS 3.720e-11, quick 27/27, full 63/63, and the payloads of the release suite are bitwise. Please build it as REL (CUDA 12.9 + g++13, star_vet, never CUDA 12.2), run the tiny gate first (pass1 must be 3.720e-11), and post the md5 and prov.

## Stars and starts
m35mms, t9l62, m35zams, t9mlt. For each, use the deck the latest production run restarted with (its own `-i` deck, with its command-line keys put into the deck) and that run's latest rst.

## The check: `docs/handover/reconfatal-1127/rstkeys.py`
```
ATHENAK=<this branch's tree> LAUNCH="<your 1-node launcher, as production>" NT=<threads> \
  python3 rstkeys.py <EXE> <RST> <DECK> <rst cycle + 3> output2/dcycle=<rst cycle + 3> output2/dt=1e30 \
  output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30 [the run's other output blocks: dt=1e30]
```
It runs 3 restarts in its working directory:
1. `a_deck`: `-r RST -i DECK`. Prints the warning lines and the list of stored keys.
2. `b_named`: `-r RST -i DECK.named`. DECK.named is DECK plus each listed key at its stored value; it is written next to DECK.
3. `c_ronly`: `-r RST` alone, as the reference.

The last line is `ROW <deck> | warnings with DECK n | with DECK.named 0 | b_named vs c_ronly: k files, maxdiff 0`. Delete the rst files the runs write. If a FATAL from 4f9931fbd stops a run, post its lines and go on to the next star.

## What to report, per star
1. The ROW line, and the md5 of the rst used.
2. The `stored` lines, each with a verdict. Compare it with the list below, which is the 3-D spherical-polar start-up of this code (the values c9805a084 puts in star_bsg / star_m35mms / star_t9l62):
```
mesh/nghost 3; mesh/layout_invariant_mode tile; time/dt_retry_grow 1.1; time/dt_retry_ngood 20;
hydro/wellbalance_dynamic true; hydro/wb_x1 true; hydro/sp_x2_periodic_image true; hydro/etotgrav true;
rad_vet/force_reference wb_arad; rad_vet/clean_init steady; rad_vet/implicit_bc_x1min flux; rad_vet/clean_maxit 200;
rad_vet/opacity table; rad_vet/vet_scatter true; rad_vet/closure sc; rad_vet/vet_gd_gl_nphi 18;
rad_vet/vet_gd_wall_interp true; rad_vet/vet_gd_wall_sym true; rad_vet/vet_gd_theta_seam_sym true;
rad_vet/time_scheme hesdirk2; problem/star_sponge_form lawson
```
   - A key with the same value is *intended*.
   - Mark *EXTRA* for a key that is not in the list, and *DIFFERS* for one whose value is not the list's. For each of those, give where its value came from, if you can tell (the production deck of an earlier code, a start-up default of the code that wrote the rst, or a hand-set key). Those are the ones the USER decides on.
3. m35zams and t9mlt: add the keys to your production decks, in one of two ways:
   - use DECK.named, which names exactly the listed keys at the stored values;
   - or run `python3 addkeys.py <deck> sp`, which adds the list above (only keys the deck does not name), then name the EXTRA ones by hand.

   Do not change a DIFFERS key's value without our answer: name it at its stored value, which keeps the run as it is. Commit the decks as you keep them now, and post the paths.
4. m35mms and t9l62: say whether the stored list equals the list above. If it does, the templates of c9805a084 are proven for them too.

Reply as `NOTE-2026-11-27zj-caltech-orion-restart-keys.md`.

-- orion (RECONFATAL, for the coordinator)
