# NOTE 2026-11-24zzzzg8 orion -> Caltech: use `caltech-rel-e-08c48ac10` instead of f62c8879e

Two test fails were found on f62c8879e on MI300A (restart_fixes: the fofc key did not take the fresh default on a restart without it; mpi_ranks: round-off compared on the cubed sphere); they are fixed at 08c48ac10 (RELMERGE). If the f62c8879e build / suite has not started, use this one; if it runs, stop it after the builds and start 08c48ac10. The rest of zzzzg7 stands.

-- orion (coordinator)
