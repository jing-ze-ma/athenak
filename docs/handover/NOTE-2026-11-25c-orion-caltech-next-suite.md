# NOTE 2026-11-25c orion -> Caltech: CUDA suite of the next line (caltech-next-79c1cabef)

Branch caltech-next-79c1cabef = the pushed release + edge / cs-layout fixes, MLT source at the stage, exact Lz_tot sum, wb_star rewrite, seam conduction, MHD (FOFC / EMF / HLLD NaN / polar FATAL), the ring-clamp check. CPU green on Orion (53/53, quick 23). Please build star_vet + the suite on idle CPUs, tiny gate 3.686e-11 (deck as for rel-f), then the full CUDA suite (53 tests, 4 lists, mpi_ranks with Slurm unset, its README count) on a free test node (02-16 first, else 02-03 / 02-09); fails at once; values against rel-f's. Then the production-restart chain for one cheap star (m35zams: +1, 25 + rst + 25 vs 50, vs rel-f 50 cycles; the NOTE lines of the new restart rule (mlt_src_time etc.)). After the HG profile; the standing order of 25b continues when nothing else is queued.

-- orion (coordinator)
