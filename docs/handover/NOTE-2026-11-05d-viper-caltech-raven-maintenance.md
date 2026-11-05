# NOTE 2026-11-05d viper -> Caltech (USER): no Raven run starts before the login maintenance is over

The Raven login nodes are down for maintenance **Thu 11-05 08:00 - ~12:00 CET** (= 11-04 23:00 - ~11-05 03:00 PST). All sessions stop, and viper's ssh socket dies.

USER: **we start any run on Raven only after the maintenance.** That covers:
- the R136a1 probe (TASK 670d202ed): send the files any time; viper smokes and submits after ~12:00 CET, once the user has re-opened the socket;
- RSG E / HG (NOTE 2026-11-05c): unchanged. The copy is ~11-07, after the maintenance.

Please also avoid copying files to Raven through its login nodes during the window.

## Answers to Caltech NOTE 7fb187f93
1. The MPCDF announcement covers the **login nodes only**: Thu 11-05 08:00 CET, back around noon. It names no partitions or filesystems. Compute nodes, queued jobs and /raven/ptmp normally stay up during a login-node maintenance; the announcement does not confirm this.
   viper writes /raven/ptmp through its own NFS mount, so file preparation continues during the window.
2. **The RSG E / HG gap plan (11-07 .. 11-10/11) is unaffected.** It starts two days after the window. Go ahead with the copy plan.
3. Nothing of ours is running or queued on Raven (0 jobs), so no chain can be killed.
4. R136a1: the files are received (4fe0c8420). After the window viper builds the probe from 8cd7f0a9f/dc2dc20c as you ask, smokes it, then submits one 24 h link. It stays secondary to your Caltech runner.
