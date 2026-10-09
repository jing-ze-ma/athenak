# NOTE from Orion: please push the AG Car ext2 opacity tables (re TASK-2026-10-09-orion-lowrho-planck)

Orion has started the low-density LTE Planck/Rosseland task (work dir /orion/ptmp/jinma/agcar_opac_1009; results
will go to branch `agcar-opac-orion-1009`). To do task step 3 (comparison with ext2's extension, delivered on the
ext2 grid) Orion needs the files from /viper/ptmp2/jinma/lbv_1008/agcar/tables_ext/:

- TABLES_EXT.md
- rosseland_ext2_* and planck_ext2_*
- the extension script (and the Ferguson 2005 GS98 / TOPS source tables it read, if small, for the step-2 overlap check)

Please push them on the data branch `agcar-opac-1009` as you offered. Until then Orion computes on its own fine grid
(log T step <= 0.05, log rho step <= 0.25) and will interpolate onto ext2 once it arrives.
