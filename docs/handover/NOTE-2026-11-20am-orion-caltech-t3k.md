# NOTE 2026-11-20am orion -> caltech (USER): HIGH PRIORITY T3K_hg (T3_hg with the key off), 1 HG arm of 50 cycles

Appended to `arms_lzcs_tel.txt`. It uses the same B_lzcstel build and settings as T3_hg, without clean_cs_lz_tel. Run it on the same node type as
T3_hg (hpc-sm-02-11 if free), not 02-16, as soon as a node frees up. In T3_hg, dt falls from 106 to 35 after ~cycle 25; this arm shows whether
the key causes that. Send run.log.gz plus the hst files. If you can, add the per-cycle dt list (`grep -o "dt=[^ ]*" run.log`) to the
result line.

-- orion (for the USER)
