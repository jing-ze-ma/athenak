# NOTE deltaai -> Delta + viper (USER): Hsml probe 3464570 got NO data (bad node); resubmitted as 3466043

- Probe job 3464570 (2 nodes gh045+gh082, ghx4-interactive, binary builds3/g2mh/Arepo md5 31b1e582) ran 19:14-19:26 CDT.
- Both runs (n576, n256) failed at launch: every task on gh082 gave `srun: Error configuring interconnect`
  (tasks 288+ / 128+); the gh045 ranks hung until the per-run timeout. Zero HSMLPROBE lines -> no answer to QUESTION 5d602691 yet.
- Cost ~3 SU. This is a node fault (same as gh011 earlier), not an Arepo problem.
- USER go: resubmitted unchanged as **3466043** with `--exclude=gh082` (2 nodes, 12 min, ~3 SU). Pending; start time unknown.
- If Delta's own diag job already answered the Hsml question, tell us and we cancel 3466043.
- job3 (g25_full) stays cancelled; resubmit only on viper go + USER go.
