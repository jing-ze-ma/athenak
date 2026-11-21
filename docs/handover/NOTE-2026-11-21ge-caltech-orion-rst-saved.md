# NOTE 2026-11-21ge Caltech -> orion: production rst copies saved (answer to NOTE fz / ga step 1)

All three production runs are stopped (USER confirmed in our session at 18:57 PST). The operator posts the per-run stop lines separately.
Before the stop, every rst present was COPIED (cp -p, never mv; size verified equal to the source) to
`/resnick/groups/carnegie_poc/jingze/prod_fresh_1121_saved/<star>/`. Nothing was deleted.

| star | file | size (B) | source mtime |
|---|---|---|---|
| bsg | bsg3d.00012.rst | 7078414947 | 18:29 |
| bsg | bsg3d.00013.rst | 7078414947 | 18:37 |
| bsg | bsg3d.00014.rst | 7078415063 | 18:46 |
| t9l62 | t9l62.00004.rst | 27988113439 | 18:20 |
| t9l62 | t9l62.00005.rst | 27988113138 | 18:31 |
| t9l62 | t9l62.00006.rst | 27988113138 | 18:41 |
| t9l62 | t9l62.00007.rst | 27988113439 | 18:51 |
| m35mms | m35mms.00001.rst | 9239544529 | 17:59 |
| m35mms | m35mms.00002.rst | 9239544830 | 18:21 |
| m35mms | m35mms.00003.rst | 9239544830 | 18:44 |

- **BSG:** rst output5 dt is 1 d, so 00012/00013/00014 are the t = 12/13/14 d dumps. Those were the oldest BSG rst still present; KEEPRST had already pruned 00000-00011.
- **Cycles and md5s:** the operator's stop lines give each rst's cycle (from the log's restart-dump lines) and md5, so you can check whether 00014 is before the cycle-9178 burst.

-- Caltech
