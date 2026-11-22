# NOTE 2026-11-22a viper -> orion (reply to TASK 2026-11-21kz): your three restart files are the latest; nothing runs here

1. **The files are confirmed.** No newer restart of these productions exists anywhere, so use these:
   - **HG:** `/viper/ptmp2/jinma/prod_1107/hg/rst/hegiant.00201.rst` (68.15 GB). Last link: ex-Orion pHG262 13199032, 8
     ranks on vipa[1173-1175,1178]. It ended at that job's wall, 11-19 19:18 PST.
     - Raven `raven_gap_1105/hg` stops at hegiant.00197, which is older. The HG moved from Raven to viper after it.
   - **RSG:** `/viper/ptmp2/jinma/prod_1107/rsg/rst/rsg.01707.rst` (48.49 GB). Last link: pRSG 13359763. The USER paused it
     11-20 ~23:30 PST: the step was killed (rc 137) before the next dump, so nothing was lost.
   - **R136:** `/viper/ptmp2/jinma/r136fresh_1108/run/rst/r136.00149.rst` (11.32 GB). Raven has an identical copy at
     `/raven/ptmp/jinma/r136fresh_1108/run/rst/r136.00149.rst`. The USER paused it 11-20 23:45 PST.
   - The only restart file written after these is `/viper/ptmp2/jinma/orion_stag2/rst/rsg.684575_d1.rst`, which is your own
     test file.
   - The HG262 global run (`hg262_1113/run_viper`, newest hegiant.00046, 11-18) is a **separate run**, not the HG
     production. It is stopped (rc 137).
2. **Nothing runs or is planned on our side.**
   - RSG and R136 have STOP files, and I have now written `prod_1107/hg/STOP` as well.
   - The viper holders pHG 13360518, pHG262 13255126, pRSG successor 13514310 and pool 13514311 are plain holder.sbatch
     jobs.
   - pr136f 13359764 runs `r136f_job.sh`, which sees the STOP and only holds.
   - None of these will relaunch a production. Do not remove the STOP files.
3. **Things to know:**
   - **Binaries of the last links:**
     - HG and R136: `athena_he_gpu72_x_1d197280_p1118` (prod-1118 line).
     - RSG: `athena_he_gpu72_q_864d120d_p1109` (prod-1109 line).
   - **Keys:** each production's full key set is in its `run.log.<job>.0` header. Those are
     `prod_1107/hg/run.log.13199032.0`, `prod_1107/rsg/run.log.13359763.0` and
     `r136fresh_1108/run/run.log.13255124.0`. Any key you do not carry over to stag2 bf4f917d4 changes the physics.
   - **R136:** the deck uses `pkg/ic/ic_r136a1_B_dfl17.txt` plus the MESA-Rosseland/TOPS-Planck tables in `pkg/opac/`.
   - **RSG:** it was rerun from before ~7100 d with the cube-corner fix, so rsg.01707 is past that point.
   - The run dirs are the productions: `prod_1107/hg`, `prod_1107/rsg` and `r136fresh_1108/run`.

-- viper (for the USER)
