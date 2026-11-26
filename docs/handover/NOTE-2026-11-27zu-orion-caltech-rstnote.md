# NOTE 2026-11-27zu orion -> Caltech: PRIORITY 1, the removed-key restart fix (your 27zr) on the merged tip caltech-next-3c681d9d1

USER decision on your 27zr (clean_retry_budget = 0 stored in every production rst FATALed):
a removed key's value that only the restart file holds (no deck, no command line names it;
also old `# Default value added at run time` lines) now prints
`### NOTE: <rad_vet/clean_retry_budget> = 0 is stored in the restart file (not named): the key was removed, the code runs true`
and the run continues with the fixed value; a value a deck or the command line names that differs
still FATALs. Commit 3c681d9d1, on top of the merged tip (release aa2646b92 + raycut, docmap,
reconfatal, auditfix (foot, retry, FEDGE), conserv (massless faces, docs), KEYS2 (8 keys removed,
47 keys renamed without clean_, old names read as aliases), REPRO (fixed user.hst order, layout
tile 16)). Branch **caltech-next-3c681d9d1** (pushed). Orion CPU: restart_keys (new arms rm_note /
rm_same / rm_fatal), restart_fixes, key_fatal, cmdline_keys PASS; full gate running.

Please, on one idle test node (never one with production), after or instead of the rest of 27zq:
1. CUDA 12.9 + g++13 build of caltech-next-3c681d9d1 (built_in_pgens + star_vet); tiny_bsg_rel5
   3.720e-11; the quick and the FULL list (27zq's request, now on this tip: it supersedes 104c68b7e).
2. Your 27zl RST-MIX arms (bsg 00009, m35mms 00011, t9l62 00012, M35ZAMS 00011), on this binary,
   **without** the workaround key `rad_vet/clean_retry_budget=true`: (a) 2 cycles straight, (b) 1
   cycle + rst + 1 cycle. Expect: rc 0, one NOTE line as above per run, no FATAL, (a) == (b)
   byte-identical final rst / hydro_w. Plus one arm per star with the key named `=true`: its final
   dump must equal (a) byte for byte (the NOTE changes nothing). Note: KEYS2's rename makes the
   production keys' new names the canonical ones; old clean_* names on the command line still work
   (one alias NOTE each).
Report FATAL / CLEAN WARNING / NOTE counts and bitwise yes/no per arm. Do not touch vet-clean-release.

-- orion (MERGE, for the coordinator)
