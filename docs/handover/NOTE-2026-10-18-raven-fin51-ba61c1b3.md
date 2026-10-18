# NOTE 2026-10-18 raven: production branch fin51-r-fix -> ba61c1b3 (fast-forward from bb360d8f)

Contents: cs-s4-1017 @ 1cd8a26e (cubed-sphere implicit VET port, stages 1-4; cs defaults as sp) + prodmerge-1018-hg (mltp-1017: opt-in
problem/mlt_closure = adaptive | pressure; ic-mesascaf-1017: tools/ic_star opacity.mesa_kap, startup re_zone hopf, 14p5M specs).
All new behaviour is OPT-IN (cs mesh, named mlt_closure); existing production decks are unaffected.
Gates vs bb360d8f (bsg/e1/hg production decks, 10 cycles):
- A100 strict-md5 (mergegate_1017 gbit 31263365): every hst/bin/rst file md5-identical, final rst data BITWISE, warnings identical.
- MI300A (Viper apudev 12265175): BITWISE.
Sites: rebuild production from fin51-r-fix ba61c1b3 when convenient (same cmake recipe); running chains need not restart.
Raven A100 binary md5 685204b0 (/raven/ptmp/jinma/mergegate_1017/bin/athena_he_a100_ba61c1b3).
He giant 14p5M production (USER go 10-18) runs on 4c494982 (= the hg half of ba61c1b3, bitwise on production decks).
