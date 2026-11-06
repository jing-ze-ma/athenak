# NOTE 2026-11-06 orion -> viper: can Orion run the 3-D Athena++ BSG reproduction? (answers TASK-2026-11-06-viper-orion-athpp-3d-repro.md)

## URGENT for Raven jobs 31975424/31975425: check the deck's line endings
`docs/handover/athenapp-bsg/linhao_input/athinput.bsg_star.rusty` on this branch has CRLF line endings on all 98 lines. Athena++ `ParseLine` strips only spaces, so `refinement = static` (no trailing comment) is read as `"static\r"`, and `multilevel` is false. **The byte-identical deck builds 768 uniform blocks with no refinement, not 2560.** Checked with this build's mesh test (`-m 2304`, e47fabb7):

| deck | MeshBlocks |
|---|---|
| CRLF as committed | 768, "Too few mesh blocks" |
| `tr -d '\r'` | 2560 |
| CRLF + `mesh/refinement=static` on the command line | 2560 |

At 2304 ranks the CRLF deck stops at startup with a fatal error, so it cannot run wrongly unnoticed there. It would run silently on the coarse grid only at <= 768 ranks. Your 160-block probe suggests your Raven copy may already be LF, but please confirm with `grep -c $'\r' athinput` and the block count in the first log.

**Open question for the USER and Caltech:** if Linhao's own copy also had CRLF on his machine, his run would have had no refinement. The CRLF may equally have come from file transfer. Other values with no trailing comment (`data_format`, `xdmf`, `dt`, `error_limit`, etc.) also keep the `\r`; I did not check each one's effect.

**This corrects item 5 and the build statement below:** the USER then allowed the build. It compiled cleanly (gcc/13, openmpi/4.1, hdf5-mpi/1.14.1; 0 errors), and `conv_profile.py` printed T 0.998, rho 0.993 (needs `module load python-waterboa/2024.06` for numpy). A speed probe is running; numbers will follow in this note.

**Answer: CPU capacity and toolchain are fine, but Orion should not replace Raven for now.** It starts in ~19 h, which is no better than Raven's range. Free disk also fits only about the first 5-15 d, not 25 d, unless the USER clears space. Nothing started, no production submitted.

**Not tested: the build and the speed.** The automatic permission check in this session blocked compiling the external clone, so there is no build and no speed probe. Items 4-5 are based on modules only, and all Orion wall times below assume Raven's per-core rate.

1. **Nodes:** p.exclusive has 98 nodes, 2 x 56 cores (112 physical, 224 threads), 512 GB per node. Max walltime 24 h per job. QoS limit is 30 nodes per user, and g25_a120 production holds 5, so **25 nodes (2800 cores) are usable**. 2560 ranks fit on 23 nodes (1 block per rank). Memory is ample: 1.1 TB total, or 2.6 TB at his 1 GB per block, against 11.8 TB on 23 nodes.
2. **Start and chains:** `sbatch --test-only` for 23 nodes, 24 h gives 2026-11-07 14:29 CET, ~19 h after this check. This is an estimate only. Restart chains are OK: g25_a120 runs as a dependency chain now. The limit is 6 running jobs per association, with 2000 submitted.
3. **Disk: the blocker for 25 d.** /orion/ptmp has a per-user quota of ~60 TB (recorded value; no quota tool works here). Usage was 57 TB on 11-05 and g25_a120 grows by ~140 GB/day, so **under ~3 TB is free**, and g25_a120 also needs a planned ~1.4 TB of snapshots. Today's total was not re-measured because `du` did not finish. ~0.9 TB per 5 d fits; 4.5 TB does not without the USER freeing space.
4. **Toolchain:** yes. gcc/13 + openmpi/4.1 + hdf5-mpi/1.14.1 modules exist (hdf5-mpi/2.0.0 reports Parallel HDF5 ON; 1.14.1 not queried). Load modules outside a pipe, or they do not take effect.
5. **Building public Athena++:** not tried, see above. e47fabb7 cloned from GitHub without trouble, and the bundle was copied to `/orion/ptmp/jinma/athpp_repro3d`.

## Cost if it ran here
At Raven's 3350 zone-cycles per core-s and 2560 ranks: ~27 h per 5 d and ~133 h (~6 links) for 25 d. That is 10 % faster than 2304 ranks on Raven at best. Orion's per-core rate is unmeasured.

## Suggestion
Keep Raven as the primary (31975424, 31975425). Use Orion only as a second copy if the USER wants one; that needs permission to build here and a decision on disk for anything past ~5-15 d. A 1-node probe can follow the build to replace the assumed rate.
