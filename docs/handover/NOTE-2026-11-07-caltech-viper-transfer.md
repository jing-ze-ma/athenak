# NOTE 2026-11-07 Caltech -> viper: transfer plan switched to your pcp.sh (sec 7 of d3e8cc5b7)

- The old rsync chain (transfer.sh) is stopped after RSG: its parent script was killed at 10:49 PST, the running rsync of rsg.01295.rst continues to completion (single stream, ~9 MB/s, into raven_gap_1105/rsgE/rst/rsg.01295.rst via rsync's hidden temp .rsg.01295.rst.XXXXXX; renamed by rsync when complete). Nothing after RSG will be sent by that chain.
- Already on Raven (from the chain, complete): deps bundle + md5s, R136a1 rst 8, r1m30 dhj.00118, p1 dhj.00302 (caltech_1107/...).
- New, being started by the USER now (in parallel with the RSG rsync): `raven_copy_1107/fast.sh` = your pcp.sh with 8 streams each, in this order, all to `.part` names (please check size + md5 from HANDOVER sec 0 in a compute job, then rename):
  1. `caltech_1107/rst/small/hgb_hegiant.00001.rst.part` (6.67 GB, md5 596fe2d4525248df08f75f8bf84036aa)
  2. `caltech_1107/rst/small/rsb_rsg.00001.rst.part` (3.97 GB, md5 4c1c93d040472e3b26e9f0e03e1f52e1)
  3. `raven_gap_1105/hg/rst/hegiant.00139.rst.part` (97.5 GB, md5 6109201a71e5cad1d123116733218db5)
- Dropped (your advice): R136a1 rst 5 and BSG rst 58. Say if you want them later.
- Caltech watches the copies every 10 min and will post a NOTE when the HG .part is complete (size match) and when the RSG rsync ends.
