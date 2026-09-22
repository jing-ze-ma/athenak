# Vertex band vs panel interior (vertex.py)

Vertex band = within 4 cells of a cube vertex in BOTH gnomonic angles; interior = the middle half of each panel in both. `res` = rms departure of the horizontal wind from its zonal mean at the same latitude and radius.

### comparison 2 (cold start)

| arm | rot | cells vert/inter | rms |v_t| vert | interior | ratio | rms residual vert | interior | ratio | grid-noise vert | interior | ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| plm_ctl | 0.250 | 384/1536 | 2.578e+05 | 3.171e+05 | 0.813 | 2.512e+05 | 3.148e+05 | 0.798 | 4.945e+04 | 5.076e+04 | 0.974 |
| ppmx | 0.224 | 384/1536 | 2.591e+05 | 3.123e+05 | 0.830 | 2.512e+05 | 3.100e+05 | 0.810 | 6.271e+04 | 6.643e+04 | 0.944 |
| wenoz | 0.215 | 384/1536 | 2.623e+05 | 3.078e+05 | 0.852 | 2.567e+05 | 3.055e+05 | 0.840 | 6.572e+04 | 7.348e+04 | 0.894 |

### comparison 1 (restart, reference)

| arm | rot | cells vert/inter | rms |v_t| vert | interior | ratio | rms residual vert | interior | ratio | grid-noise vert | interior | ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| prodbin | 284.269 | 384/1536 | 4.799e+05 | 4.862e+05 | 0.987 | 4.578e+05 | 4.815e+05 | 0.951 | 3.938e+04 | 7.082e+04 | 0.556 |
| sph | 283.736 | 384/1536 | 4.373e+05 | 4.695e+05 | 0.931 | 4.200e+05 | 4.632e+05 | 0.907 | 5.398e+04 | 5.592e+04 | 0.965 |
| sphbeam | 283.648 | 384/1536 | 5.010e+05 | 4.711e+05 | 1.064 | 4.824e+05 | 4.646e+05 | 1.038 | 3.645e+04 | 5.812e+04 | 0.627 |
| noang | 283.735 | 384/1536 | 4.373e+05 | 4.711e+05 | 0.928 | 4.201e+05 | 4.641e+05 | 0.905 | 5.390e+04 | 5.308e+04 | 1.016 |

### comparison 2: vertex/interior residual ratio vs radius

| r/Rp | plm_ctl | ppmx | wenoz |
|---|---|---|---|
| 1.004 | 0.60 | 0.76 | 0.92 |
| 1.088 | 0.77 | 0.89 | 0.98 |
| 1.151 | 0.98 | 0.94 | 0.80 |
| 1.210 | 0.84 | 0.79 | 0.67 |
| 1.275 | 0.76 | 0.84 | 0.66 |
| 1.352 | 0.72 | 0.70 | 0.80 |
| 1.445 | 1.11 | 1.21 | 1.17 |
| 1.557 | 0.40 | 0.39 | 0.47 |
| 1.690 | 0.50 | 0.48 | 0.52 |
| 1.848 | 0.82 | 0.79 | 0.87 |
| 2.038 | 1.14 | 1.10 | 1.22 |
