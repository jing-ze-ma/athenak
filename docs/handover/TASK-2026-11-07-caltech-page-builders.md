# TASK 2026-11-07 viper -> Caltech (USER: "do you have all the tools to update rsg hg bsg r136 urls? if not ask caltech to send you")

Viper now runs these productions (Raven/viper) and must keep their pages current. Viper can publish to the existing artifacts (same account) but does NOT have your CURRENT page builders (the branch only has October versions: rsgE-page-deltaai, hegiant-page-1007, bsg-hrdet-1009/page, bsgprod-1012/pagework).

Pages to take over (please confirm each URL is the live one):
| run | page | URL |
|---|---|---|
| RSG E | "RSG E Caltech" | https://claude.ai/artifact/LBASg2iK2ArcLY1V2yBCAR |
| HG | "He Giant Global Caltech" | https://claude.ai/artifact/LNtjtaL34oUhgcCRNp5V2M |
| R136a1 B | "R136a1 deep envelope run (B)" | https://claude.ai/artifact/J8TiLW4e42UMvjWCqQyvRh |
| BSG | "BSG Linhao-exact (production)" (prod_lhxo, page v24) | https://claude.ai/artifact/GyzAhm3DR7F81DyaBfvL2n |
(+ if a W121 r1m30 page exists, e.g. "Hot Jupiter Round 1" https://claude.ai/artifact/A5omHrtMwF5roVtjNjDqq7 — say which.)

Please push to this branch, under `docs/handover/caltech-to-raven-1107/pages/<run>/` (one dir per page):
1. The exact current builder scripts (python + shell), unchanged, and a README with: the commands in order, the inputs each step reads (hst files + which columns, which dump types/cadence, derived caches), what is computed (mode analyses, slices, movies, budgets), the output html (+ any asset files and how they are published: inline data URIs vs Artifact assets), and the page's current version number.
2. Any CACHE / derived-data files the page needs to keep its history continuous from the Caltech part (e.g. time series derived from dumps that are not in the hst; per-snapshot extracts). If larger than ~50 MB, put them in `$J/raven_copy_1107/pages_cache.tar` + md5 and say so; the USER's fastpush route can send it to /raven/ptmp/jinma/caltech_1107/pages/.
3. Machine-specific bits to repoint (absolute $J paths, module names, matplotlib backend, fonts).
4. Any page conventions the USER asked for (layout, units, which plots first, wording rules) that are not in the scripts.
Viper will then rebuild each page from the Raven/viper outputs (appending to your history) and republish to the same URL; first rebuild without new data, to check the page comes out identical before new data is added.
