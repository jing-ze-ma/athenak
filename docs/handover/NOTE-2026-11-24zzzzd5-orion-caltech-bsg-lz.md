# NOTE 2026-11-24zzzzd5 orion -> Caltech: bsg: does the merge fall-back cost L_z? (read-only, no run)

The merge fall-back (nfb) sets a group's rotation to zero: M, E and linear momentum stay, the group's spin about its centroid is dropped. c1 bsg: 0 to c6000, 101,260 at c7000, 575,785 at c9000. From existing hst and logs only:
1. c1 bsg (and the T2 bsg, and rc3 bsg if it has Lz_tot): Lz_tot(t) and its relative drift per 500 cycles; the same for the angular-momentum columns the hst has (other components if any). Does the drift per cycle change at c6000-9000, where nfb starts? Compare with bsg at cycles 0-5000.
2. Scale: |L_z| of the whole star, its initial value, and the L of the gas in the merge region (if a column or dump allows).
3. MLATDIAG lines: the L column drift over the same window, and the "chk" columns if printed.
4. Where are the fall-back groups: is there any per-cell print? If not, say so; the radius range of the merged groups with Erot > half Eth can be read from the newest bsg dump (bin) if one exists: count groups (2x2 cells, r < 2.43e12 f = 2, 3.83e12 f = 1, 5.57e12 f = 2) with 0.5 omega.L > 0.5 e_int, their rho and e_int range. CPU, low priority, not on a production node.
5. Don't change anything. bsg goes on.

-- orion (coordinator)
