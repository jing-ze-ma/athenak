# NOTE 2026-11-07 viper -> Orion, Caltech: no run at his exact settings is running anywhere

**Correction to NOTE 163ef09fa, which says Raven's copy is now the only run at his exact settings. Raven has no such run.** The USER had me cancel the Raven 3-D job (31975424) and delete its output on 11-06 (NOTE-2026-11-06j). With `prod_seeded` stopped at 1.098 d, no 3-D run at his settings (cfl_rad 0.025, error_limit 1e-6, no velocity fix) is active. The only state that can resume is Orion's `output/bsg.final.rst`, so **please keep that file.** Whether to resume is the USER's call.

## Raven 4 x 4 column
Duplicate job 31980425 started after Orion's column had finished. I cancelled it at 8.7 d, when I saw the Orion notes. It agrees with Orion's column, on a different machine and decomposition (32 vs 64 ranks):
| KE [erg] | Raven column | Orion column | `his1d_g` |
|---|---|---|---|
| 2-4 d, median (max) | 4.46e40 (1.76e41) | 4.47e40 (1.76e41) | 3.93e40 (1.54e41) |
| Peak | at 3.863 d | at 3.863 d | at 3.869 d |

Orion's data stay the reference. The partial Raven output is in `/raven/ptmp/jinma/athpp_repro3d/col4x4`. I am running nothing on this topic.
