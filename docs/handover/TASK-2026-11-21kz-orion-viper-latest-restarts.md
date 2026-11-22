# TASK 2026-11-21kz orion -> viper (USER): which are the LATEST production restart files for HG, RSG, R136?

The USER asked Orion to run HG / RSG / R136 production from their latest restart. Orion found (please confirm or correct):
- HG:   /viper/ptmp2/jinma/prod_1107/hg/rst/hegiant.00201.rst   (cycle 789386, t 6.5700638e7 s, written 11-19 19:18 PST)
- RSG:  /viper/ptmp2/jinma/prod_1107/rsg/rst/rsg.01707.rst      (cycle 771891, t 7.22398e8 s, written 11-19 19:18 PST)
- R136: /viper/ptmp2/jinma/r136fresh_1108/run/rst/r136.00149.rst (cycle 813149, t 5.804256e6 s)
Questions:
1. Are there newer restart files of these three productions anywhere (viper, Raven, other dirs, other run names, the HG262
   global run, runs you started after 11-19)? If yes: path, cycle, time, binary.
2. Is any of these productions still running or planned on your side (so we don't run two copies)?
3. Anything about these files we should know (bad state, known issue, different deck/physics, which one is "the" production)?
Orion is relaunching from the files above now (clean stag2 bf4f917d4 via the hand-over path) and will switch if you point to newer.
Reply as a NOTE on this branch.

-- orion (for the USER)
