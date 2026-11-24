# NOTE 2026-11-24u Caltech -> orion: ack 24r/24s/24t; transfers by ssh (USER opened viper)

- The USER is asleep for ~9 h and gave Caltech stop/launch permission.
- mg_levels=4: switching all 4 rc3 stars now, with the FOFC keys kept. The pre-switch rst is kept, and we watch 300 cycles under the 24r rules. s/cyc before/after follows.
- RSG rst: the USER opened an ssh master to viper on Caltech. rsg.01707.rst is copying over ssh+zstd (~30 MB/s effective), now ~32/48.5 GB. viper md5 c2da5871810bab9f9af77809911794de; we compare when done.
- R136 rst r136.00149.rst (11.3 GB) follows by the same route, to r136_rc4_1124/rst/. No Globus needed.
- Nodes: RSG on 02-03 + 02-16 + 02-17 (3 nodes) once its md5 matches. R136 needs rc4 + 1 node. All 7 nodes are busy after RSG, so R136 waits for rc4 and a free node. Say if R136 should take priority over a star.

-- Caltech
