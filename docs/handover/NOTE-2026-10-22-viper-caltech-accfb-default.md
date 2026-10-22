# NOTE 2026-10-22 viper -> caltech (USER): time2_accept_fallback default-on = accfb-default-1022 1976bac68; fin51-r-fix fast-forwarded

- accfb-default-1022 = 1976bac68 = ff1f52bc + <rad_m1>/time2_accept_fallback DEFAULT true under hesdirk2 (key false turns it off).
  A fresh -i run cannot name the key on the command line unless it is in the deck (absent-key override FATALs).
- fork/fin51-r-fix fast-forwarded ff1f52bc -> 1976bac68.
- MI300A gates (viper apudev, 2 ranks):
  - RSG E det, HG global det (cs), hg145 (sp wedge): fallback 0 -> bitwise vs ff1f52bc.
  - RSG E C2==C, R==C; key-off decks bitwise vs ff1f52bc (RSG E, accretor D_B1v6).
  - Accretor D_B1v6 static to 4e-4: vimp fallback fires; default-on 63 cycles / 6 retries vs 64 / 7 (cycle-57 retry removed); C2==C, R==C.
  - Accretor stream (st128, extrap2 + harmonic, t 0.004 -> 0.0055): on == off == ff1f52bc bitwise (fallback only in a non-admissible stage).
  - No rad_m1 MHD deck in the gate set (not covered).
- Production smoke viper 12328922 (run12 deck, static + stream restart) clean, accepted-with-fallback stages=1.
- Builds: MI300A accretor md5 c75cbac9, MI300A he_star_m1 md5 ef3a984f, A100 accretor nofma md5 9b194823.
- viper run12 (apu, pending) now runs 1976bac68. Raven run12 + HG wedge switch pending.
- Next: posfix for am_vector (TASK 92c087ad) taken.
