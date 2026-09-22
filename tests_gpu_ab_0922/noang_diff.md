# noang vs sph (noang_diff.py)

hst rows compared: 45, times identical: False

| column | max |rel diff| | first row with rel diff > 1e-10 (t, rot since restart) |
|---|---|---|
| time | 2.311e-06 | 8.653160e+07 (0.4315) |
| dt | 6.746e-03 | 8.651630e+07 (0.3813) |
| mass | 0.000e+00 | never |
| 1-mom | 6.737e-03 | 8.641260e+07 (0.0413) |
| 2-mom | 4.039e-04 | 8.648280e+07 (0.2715) |
| 3-mom | 4.046e-03 | 8.643700e+07 (0.1213) |
| tot-E | 0.000e+00 | never |
| 1-KE | 2.121e-03 | 8.640960e+07 (0.0315) |
| 2-KE | 1.849e-03 | 8.644310e+07 (0.1413) |
| 3-KE | 3.266e-04 | 8.640960e+07 (0.0315) |
| 1-ME | 1.645e-04 | 8.643400e+07 (0.1115) |
| 2-ME | 1.005e-04 | 8.645840e+07 (0.1915) |
| 3-ME | 3.134e-04 | 8.643090e+07 (0.1013) |

dt: identical in every row: False

last dumps: sph t=8.654230e+07 cycle 4440810, noang t=8.654210e+07 cycle 4440789 (NOT the same instant, so a cellwise diff mixes the physics difference with 200 s of evolution)

| variable | rms rel diff | max rel diff | location of max (lat, lon, i) |
|---|---|---|---|
| dens | 4.894e-05 | 1.129e-02 | +9, -156, 17 |
| velx | 9.924e-02 | 6.531e+00 | -7, -176, 55 |
| vely | 1.786e-01 | 7.810e+00 | -7, -173, 52 |
| velz | 6.006e-02 | 3.866e+00 | -27, -176, 53 |
| eint | 3.675e-05 | 7.242e-03 | +9, -156, 18 |
| bcc1 | 4.816e-02 | 6.629e+00 | +12, -162, 25 |
| bcc2 | 3.705e-02 | 4.536e+00 | +16, -150, 17 |
| bcc3 | 4.142e-02 | 5.099e+00 | +16, -153, 40 |
