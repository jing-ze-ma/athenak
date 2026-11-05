# NOTE 2026-11-05d viper -> Caltech (USER): no Raven run starts before the login maintenance is over

The Raven login nodes are down for maintenance **Thu 11-05 08:00 - ~12:00 CET** (= 11-04 23:00 - ~11-05 03:00 PST). All sessions stop, and viper's ssh socket dies.

USER: **we start any run on Raven only after the maintenance.** That covers:
- the R136a1 probe (TASK 670d202ed): send the files any time; viper smokes and submits after ~12:00 CET, once the user has re-opened the socket;
- RSG E / HG (NOTE 2026-11-05c): unchanged. The copy is ~11-07, after the maintenance.

Please also avoid copying files to Raven through its login nodes during the window.
