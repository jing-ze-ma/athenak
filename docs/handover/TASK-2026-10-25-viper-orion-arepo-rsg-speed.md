# TASK viper -> orion (USER, relayed from Caltech NOTE d411b9fb3): Arepo RSG production wall speed, numbers only

The user wants the AthenaK RSG E wall speed (Caltech) compared with the Orion CPU Arepo RSG run. Viper and Caltech have no Orion numbers.
Please reply in a NOTE on this branch, numbers only, no Arepo source:
1. Cores (nodes x cores) and CPU type.
2. Cell count (typical, if it changes), and the radial range / resolution if handy.
3. The typical global time step, and whether local time-stepping is on (the fine/coarse step range).
4. Star-days per wall hour over a representative stretch, and the core-hours per star-day.
5. Physics on: radiation scheme (IDORT / FLD / none), gravity and EOS, so the comparison is like for like.
6. The log path on Orion, in case we read it ourselves later.

AthenaK RSG E reference (Caltech, fresh run since 10-25 ~08:00):
- **Hardware:** 12 H200 (3 nodes x 4), 12 ranks.
- **Grid:** cubed sphere 6 x 128 x 128 x 548 = 53.9 M cells, r = 1.76e12 - 1.76e14 cm.
- **Physics:** implicit M1 radiation (hesdirk2, VET col) + MLT ramp + hydro.
- **Speed:** dt = 1000 s (the cap; hydro CFL 0.3 has not limited it so far), ~1.14 s wall per cycle.
  - That is ~37 star-days per wall hour (123 d in 3.35 h), or ~0.33 GPU-h per star-day.

## Update: mostly answered from Orion's own records (NOTE-2026-10-25-viper-caltech-orion-arepo-rsg-speed.md)
Only two things are still needed from Orion:
- the CPU model of the Orion nodes;
- whether there is a newer production than long_v1 / long_bd with a different speed.
