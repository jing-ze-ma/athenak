# NOTE 2026-11-24p orion -> Caltech: t9l62 FOFC ON (USER); RSG transfer route

## 1. t9l62: turn FOFC on (USER, after your WD0 / WD2)
- At the next t9l62 rst of the fresh rc3 run, restart with the same deck / args plus `hydro/fofc=true hydro/fofc_retest=true`.
- Check the log: `<hydro> ... fofc_retest 1`. The rst stores fofc = false, so name both keys.
- Don't trust "FOFC first-order cell-stages": that counter is broken (reads 0 while FOFC acts). Use the end-of-run "energy created x vol" line and the floor_audit pay lines instead. A counter fix is coming.
- Post s/cyc before / after the switch (FOFC's cost), CW, dt.
- The other stars stay as they are until the USER decides (we're measuring the cost).

## 2. RSG rst transfer (48.5 GB on viper)
Orion has no ssh route to Caltech, and viper's /ptmp2 isn't mounted on Orion. So the clean route is **Globus**, started once by the USER (one 2FA login), unattended after that:
- source collection: MPCDF's Globus endpoint for viper's /ptmp2 (search "MPCDF" in the Globus web app), path `/viper/ptmp2/jinma/prod_1107/rsg/rst/rsg.01707.rst`;
- destination collection: Caltech HPC's Globus collection, path `/resnick/groups/carnegie_poc/jingze/rsg_rc3_1124/rst/`.
If either side has no Globus collection, tell us and the USER decides on the gate login.

-- orion (coordinator)
