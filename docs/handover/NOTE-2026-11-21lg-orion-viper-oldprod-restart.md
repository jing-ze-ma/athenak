# NOTE 2026-11-21lg orion -> viper (USER): production reverted to the OLD binaries; HG stopped

USER (11-21 ~19:45 PST) stopped the clean-scheme prod1121 runs and reverted to the old production binaries from the pre-switch rst
(your NOTE 2026-11-22b jumps: confirmed as a numerical artefact of the clean scheme's thin floor-gas transport; being fixed).
Running since 18:29 PST (Orion clock), in NEW dirs (your prod_1107 / r136fresh_1108 dirs and STOP files untouched; read-only):
- RSG: viper pRSG 13359763, /viper/ptmp2/jinma/orion_stag2/oldprod1121/rsg, athena_he_gpu72_q_864d120d_p1109 (md5 9acb19df),
  from rsg.01707 (hardlink), keys = cmd.13359763 exactly; driver oprun.sh (DSTOP / STOP_NOW), relaunches on successor 13514310.
  Please make sure 13514310's holder script does not start a run of its own (it should see STOP + HOLDER_ON in prod_1107 and hold).
- R136: Raven pool 32638311 (A100), /raven/ptmp/jinma/stag2/oldprod1121/r136, Raven CUDA build of 1d1972805 (md5 b77d9eb1),
  from r136.00149, keys = cmd.13255124 (table paths only changed); successor 32747999.
- HG: STOPPED (USER); last rst /raven/ptmp/jinma/stag2/prod1121/hg/rst/hegiant.00204.rst (clean scheme, cycle 791620).
Your pages (HG/RSG/R136) can follow the oldprod1121 dirs for the new segments; mark the clean segment (16:10-19:45 PST) as such.

-- orion (for the USER)
