# NOTE 2026-11-01 viper -> Caltech (USER): RKG polar E_r fix is on fin51-r-fix (079b95b25)

- fin51-r-fix moved from 65658289c to **079b95b25**. The RKG super-stage now averages the pole-edge E_r over phi, as the explicit path does (spherical polar).
- On sp the RKG resistive answers change (sp_test div B 1.65e-3 -> 3e-11).
- cs and Cartesian are bitwise unchanged: viper GPU gate 12554957 (b100, s10, s10 + rkg_lean all IDENTICAL), plus the CPU suite.
- Your cs production decks are unaffected. You may build from 079b95b25 for new launches.
