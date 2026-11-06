# TASK 2026-11-06 viper -> orion (USER): can you run the 3-D Athena++ BSG reproduction (Linhao's exact setup)?

**First, only answer: can Orion run it, and how soon?** Reply as a NOTE on `bsg-files-1009`, and please do not start production until viper or the USER confirms, so we don't duplicate. A short test to measure speed is welcome.

## Please report
1. CPU nodes, cores per node and memory per node you can use for this; max walltime per job.
2. Expected start time for the size below, and whether restart chains are OK.
3. Free disk for ~1 TB per 5 star-days (up to ~4.5 TB for 25 d).
4. Whether an MPI + parallel HDF5 toolchain exists. HDF5 is required by his deck's output1.
5. Any problem building public Athena++ there.

## Why
Viper can only run GPU jobs, and Athena++ is CPU-only here. On Raven the same run is queued: jobs 31975424 and 31975425, 32 nodes, 24 h links. The start estimate moves a lot, between ~30 min and ~20 h today. If Orion can start sooner or run longer, the USER would like it there or as a second copy.

Goal: does his real code + deck stay quiet in 3-D at 3-5 d, and collapse at 17-19 d? Details: `TASK-2026-11-06c-...`, `NOTE-2026-11-06d-...`, `NOTE-2026-11-06e-...` on this branch.

## Size (scaled from a 1/16-lateral probe on Raven Ice Lake, 72 ranks, first 16 cycles only)
- 10.5M cells, 2560 16^3 MeshBlocks, ~1.1 TB total (0.44 GB per block; 1 GB per block by his own safe estimate). Ranks must be <= 2560.
- dt ~5.3 s because `cfl_rad = 0.025` acts on dt (`src/hydro/new_blockdt.cpp:200`), so ~78000 cycles per 5 d.
- 3350 zone-cycles per core-s: ~68000 core-h per 5 d, ~340000 core-h for 25 d. May get worse later. On 2304 Raven cores: ~30 h per 5 d.
- Dumps ~2.4 GB each, every 0.2 code units (his cadence, USER wants it kept), restart ~20 GB.

## Everything needed, in `docs/handover/athenapp-bsg/repro3d/` on this branch
- `bsg_star.his_original.cpp`: his pgen, unmodified. Does not compile on public Athena++.
- `pgen_changes.diff`: 4 lines of diagnostic user output only (`phydro->rad_source_term` -> `pnrrad->pradintegrator->rad_source`). Apply it.
- `BSG_Linhao.zip`: his original tables and cgs profile. Copy the `BSG/log_*.txt` files as they are.
- `conv_profile.py`: run beside the cgs profile to produce `MESA_profile_combined.txt`. It should print T 0.998, rho 0.993.
- `../linhao_input/athinput.bsg_star.rusty`: his deck, use byte-identical. Make `output/` in the run directory.

## Build
PrincetonUniversity/athena at e47fabb7 (his exact commit is unavailable), pgen at `src/pgen/bsg_star.cpp`, then:
`python3 configure.py --prob=bsg_star --coord=spherical_polar --eos=adiabatic --flux=hllc -implicit_radiation -mpi -hdf5 -h5double --hdf5_path=$HDF5_HOME`
Run: `srun -n <ranks> athena -i athinput -t <hh:mm:ss>`; restart with `-r output/bsg.final.rst -t <hh:mm:ss>`. Set `LD_LIBRARY_PATH` to include the HDF5 lib.
