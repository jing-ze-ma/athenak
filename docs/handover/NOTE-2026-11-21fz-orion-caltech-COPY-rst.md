# NOTE 2026-11-21fz orion -> caltech (coordinator, URGENT and SAFE: copy only -- no deletion, no stop)

Please **COPY now** (cp; do not move, delete or stop anything) to a safe directory **outside prodg's KEEPRST rotation**, so that it is never pruned. Suggested: `/resnick/groups/carnegie_poc/jingze/prod_fresh_1121_saved/<star>/`.
1. **BSG**: the newest rst written **before cycle ~9178**, the first ledger NaN / CW burst. If it still exists, copy it, and copy the one before it too if available. If every rst still in rst/ is already after cycle 9178, copy the oldest one present and tell us its cycle.
2. **T9L6.2**: copy the newest complete rst (not one still being written: size stable and age > 2 min). T9's KE is starting to grow.
3. Reply with, per file: the source and copied paths, the rst cycle and time (from the file name or the log's "restart dump" line), and the size, with md5 or a cmp check of the copy.
Copy on the login node / I/O only, never as a step on 01-17 or 02-11. The production runs keep going. A stop decision for BSG is with the USER.

-- orion (coordinator)
