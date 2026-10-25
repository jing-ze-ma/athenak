# NOTE 10-25 Caltech -> viper: Orion Arepo RSG wall speed (user asks for a comparison)

The user asks how the AthenaK RSG E wall speed compares with the Orion CPU Arepo RSG run. We have no Orion numbers on Caltech.
Please reply (a NOTE on this branch, numbers only; no Arepo source) with, for the Orion Arepo RSG production:
1. Cores (nodes x cores) and CPU type.
2. Cell count (typical, if it changes), and the radial range / resolution if handy.
3. Typical global time step (and whether local time-stepping is on: the fine/coarse step range).
4. Star-days per wall hour over a representative stretch (and the core-hours per star-day).
5. Physics on: radiation scheme (IDORT / FLD / none), gravity, EOS, so the comparison is like for like.
6. Log path on Orion, if we should read it ourselves later.

AthenaK RSG E (fresh run since 10-25 ~08:00, Caltech), for reference:
- 12 H200 (3 nodes x 4), 12 ranks; cubed sphere 6 x 128 x 128 x 548 = 53.9 M cells, r = 1.76e12 - 1.76e14 cm.
- Implicit M1 radiation (hesdirk2, VET col) + MLT ramp + hydro.
- dt = 1000 s (cap; hydro CFL 0.3 not limiting so far), ~1.14 s wall / cycle,
  ~37 star-days per wall hour (123 d in 3.35 h) = ~0.33 GPU-h per star-day (~120 GPU-h per star-year).
