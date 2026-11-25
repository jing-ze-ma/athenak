# NOTE 2026-11-24zzzzd8 orion -> Caltech: cs arm L: where is the missing KE? (read-only, no run)

L passes the pattern rule (m4 x55 smaller) and L_top, but fails the KE rule I fixed (KE_int -15.4 % against 1 %); I don't change the rule. KE_int of Bp contains the pattern's own kinetic energy, so the question is whether the 15 % is the pattern. From Bp's and L's c600 hydro_w (and, if you kept them, c150 / 300 / 450; KE from hst too):
1. KE_int by component (radial, lateral), by radius in the imprint bins, and by cube-symmetric part (the mean over the 48 symmetry images, or at least the 6-panel mean field) against the rest: (a) KE(Bp, symmetric), (b) KE(L, symmetric), (c) KE of the rest in both, (d) KE(Bp) - KE(L) per bin.
2. KE_int(t) of both from hst (every 100 cycles): is the 15 % gap constant, growing (as t^4 for a constant force: v ~ t, KE ~ t^2 or pattern ~ t^2 -> KE ~ t^4), or at the surface?
3. L_bot, L_mid, M_tot, Etot, Lz_tot of both, relative differences at matched times, and |v| max.
4. Photosphere 800-850 Rsun KE: m4 there fell by 24 %.
Low priority, CPU. If a symmetry image tool does not exist, do the 6-panel mean and say so.

-- orion (coordinator)
