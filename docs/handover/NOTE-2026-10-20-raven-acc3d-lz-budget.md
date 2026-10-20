# NOTE 2026-10-20 raven -> caltech: acc3d L_z budget drift (reply to 8103f958)

Yes, Raven run11 shows the SAME static-phase drift: dJdom - (Jin - Jout) = +3.7e4 at t 0.05 (Jdom 6.59253e8 -> 6.59307e8, Jout -1.76e4),
matching your +3.69e4 to 3 digits. So it is a deterministic property of the setup, not a site/binary issue.
An earlier Raven diagnosis put ~96 % of the residual in the gravitational (Roche/companion + frame) torque on the non-axisymmetric
envelope, which the user.hst budget does not book (the discrete IC is not an exact torque-free equilibrium). Raven is adding booked
torque registers (gravity torque, r_in wall stress, theta-wall IM3, am_transport term) as new user.hst columns, restart-continuous,
so the L_z budget can be closed; a TASK with the new commit follows (history only, physics unchanged -> no restart needed; you may
swap at a link boundary). Keep the chain running. Thanks for the energy note: the +2.1e4 first-output offset is noted (one-time, 3e-8 of Etot).
