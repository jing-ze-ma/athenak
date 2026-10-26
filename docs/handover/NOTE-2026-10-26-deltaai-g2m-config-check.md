# NOTE 2026-10-26 DeltaAI -> viper + Delta: cfg2/g2m Config.sh checked against viper 9c28628b item 3; job3 still HELD

What I compared: DeltaAI's `cfg2/g2m/Config.sh`, the file Delta built from, against `builds/g25g/Config.sh`, the 10-24 recipe g25g built here on 10-25. The comparison is by option name, ignoring comments and blank lines.

| option | Raven g2m (viper's list) | DeltaAI g2m |
|---|---|---|
| FORCETREE_NODELOCAL | added | added |
| FORCETREE_NODELOCAL_GPU | added | added |
| IDORT_GPU_VERIFY | removed | removed |
| OUTPUT_ACCELERATION | added | **missing** |

- **Result:** nothing else differs. The only mismatch with Raven is that DeltaAI's g2m lacks OUTPUT_ACCELERATION. It is an output option, so it should not change the density loop or MaxDelaunayRadius. On config alone, I see no suspect.
- **Not covered:** values with a parameter (for example the VCELL_SEED_FROM_MAXDELAUNAY factor) count as equal only when the whole token matches. Since none showed up in the diff, they are identical.
- **job3 3462816:** still HELD (JobHeldUser), 0 SU. Release only on viper's go.
