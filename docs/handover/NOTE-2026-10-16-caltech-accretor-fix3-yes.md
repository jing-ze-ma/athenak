# NOTE 2026-10-16 caltech -> raven: accretor donor light, fix 3 = YES (USER), reply to 3bd0da2f

USER decision (10-16): please run your proposed **fix 3 test**, the vet_col lateral-beam check. vet_col builds the
Eddington tensor from radial columns, so light arriving from the side (azimuthal donor light) may not be represented.
The question is whether part of the 2-D ~30x over-irradiation comes from the closure rather than from the band
geometry, so that we know whether the donor-light scale is covering a solver limitation.

- Test only, short. Design is your choice. Suggestion: a lateral beam (or the current donor_irr deck) compared with a
  reference closure where available (M1 chi(f) / an angular SC check), with the azimuthal F/cE and E measured at r 9-10.
- Report back by NOTE: does vet_col mis-handle lateral flux, and by how much.
- Fixes 1 + 2 (no donor light from the back-sink sector; a diluted donor source) are being implemented at Caltech
  together with the user's further choices: the 2-D donor_irr scale, von Zeipel gravity darkening of the lobe
  surface, and a back-sink sector holding the stream's own light, not donor light. They arrive as a v4 TASK.
- Fix 4 (a force/sponge guard for the atmosphere) is not decided.
