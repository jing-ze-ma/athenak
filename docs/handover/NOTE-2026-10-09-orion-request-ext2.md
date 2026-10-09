# NOTE 2026-10-09 (orion -> viper): please push the ext2 opacity tables

Orion is executing TASK-2026-10-09-orion-lowrho-planck (LTE Planck/Rosseland means at log T 3.45-4.3,
log rho -21..-12, X 0.36 Y 0.62 Z 0.02 GS98). The ext2 tables exist only on viper
(`/viper/ptmp2/jinma/lbv_1008/agcar/tables_ext/`). Please commit to a data branch **`agcar-opac-1009`**
and push to origin:

- `TABLES_EXT.md`
- `rosseland_ext2_*`, `planck_ext2_*` (all files)
- the extension script that built them
- if small: the Ferguson 2005 GS98 Rosseland + Planck tables for X 0.35/0.4, Z 0.02 that ext2 used
  (orion has the TOPS x0.36 z0.02 tables but no Ferguson tables)

Orion's results go to branch `agcar-opac-orion-1009` (this branch). They are on a grid fine enough to be
interpolated onto ext2 (log T step <= 0.05, log rho step <= 0.25); a ready-to-run comparison script is
included there if ext2 is not pushed before orion finishes.
