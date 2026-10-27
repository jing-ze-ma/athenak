# NOTE viper -> caltech (USER): vet_gd theta wall, FINAL viper results and the recommended run12 option (not a GO yet)

Follows de4be1b63 (the bug) and 01fee5047 (the correction: closed hse walls, so mirror). All viper work on this is finished.

## Branch
`fork/wallbug-1027` = **917a084cf**, on 6fac5c7ff. New key: `<rad_m1>/vet_gd_theta_wall = periodic | mirror | vacuum | extrap`.
- The key is read only when named. The default `periodic` is bitwise to 6fac5c7ff on the accretor deck (hydro_w, m1, m1_face, m1_vet) and costs +0.7 % (noise).
- Every option is GPU run-to-run bitwise.
- Refused combinations: cubed sphere, `vet_gd_halo_ipc`, and extrap with more than one theta MeshBlock. run12 has one theta block.

## Accretor run12 deck (viper, X6 restart, cycles 682 -> 882, production keys)
Changes in the outer thin region (i 300-420) of the wall rows, against periodic:

| option | E change | BE steps / stage fallbacks | cost |
|---|---|---|---|
| mirror | about 5 % (median) | as periodic | about +6 % |
| vacuum | +959 % | 1 / 0 | +4 % |
| extrap | -95 % | 6 / 5 | +14..34 % |

- The differences saturate within about 50 cycles.
- vacuum and extrap are for walls open to empty space; they do not apply to run12.

## Recommendation from viper (USER decision pending; do not switch yet)
- `vet_gd_theta_wall = mirror` together with `vet_gd_gl_nmu = 11`, which makes the direction set mirror-symmetric.
- With the production GL10x12 staggered set, mirror is only approximate: the error near the phi walls grows with resolution (p90 0.10 -> 0.43 on the CPU test). With nmu odd or `vet_gd_gl_stagger = false` it is exact (wall layer within 0.1 % of the mirrored domain).
- If the user says GO, a TASK will follow. It will give the binary hash, the two key lines, a cost check (nmu 11 vs 10) and the restart point. Until then, keep run12 as it is.
