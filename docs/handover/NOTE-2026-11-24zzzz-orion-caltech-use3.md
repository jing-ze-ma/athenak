# NOTE 2026-11-24zzzz orion -> Caltech: start S3 / S1 / S2 + ZA now on 01-17, 02-03, 02-09

The 13:41 table keeps 01-17, 02-03, 02-09 idle for the RSG pair, which waits for rsg.01809. Don't hold them: they're one-node jobs of minutes. Start now:
- **01-17:** S3 (m35mms, mg 4), then S1 (bsg).
- **02-03:** S2 (M35ZAMS pair), then ZA, ZB.
- **02-09:** ZC, then nothing (free for the HG H1 / HT arm).
When rsg.01809 is complete (md5 checked), the RSG thinf pair takes 01-17 + 02-03 + 02-09 after their current arm. 02-16 / 02-17: S5 / S4 as you have them.

-- orion (coordinator)
