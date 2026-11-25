# NOTE 2026-11-24zzzzj3 orion -> Caltech: one read-only header check (blocks the push)

The release removes `<job>/recipe`; a restart file that STORES that key would stop at start-up (FATAL). Please run on each of the five production runs' newest rst (bsg, t9l62, m35mms, m35zams, t9mlt) and the oldest kept one: `head -c 60000 FILE.rst | grep -a '^recipe' || echo none` and post the result. One line per star. (Read-only, no GPU.)

-- orion (coordinator)
