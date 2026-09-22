# T(p) by region (tp.py)

Horizontal mean of T interpolated to each pressure level; day = |lon_ss|<30 deg, terminator = 75-105 deg, night = |lon_ss|>150 deg. Each arm at ITS OWN last dump (times differ, see cfl_tables.md).

### day: T [K]

| p [bar] | prodbin | sph | sphbeam | noang | max spread |
|---|---|---|---|---|---|
| 1e-06 | 4321 | 4246 | 4279 | 4246 | 76 |
| 1e-05 | 4156 | 4059 | 4088 | 4059 | 97 |
| 1e-04 | 3826 | 3695 | 3706 | 3695 | 131 |
| 1e-03 | 3476 | 3397 | 3406 | 3397 | 79 |
| 1e-02 | 2953 | 2649 | 2656 | 2649 | 305 |
| 1e-01 | 2763 | 2734 | 2737 | 2734 | 29 |
| 1e+00 | 3343 | 3341 | 3341 | 3341 | 2 |

### term: T [K]

| p [bar] | prodbin | sph | sphbeam | noang | max spread |
|---|---|---|---|---|---|
| 1e-06 | 3102 | 2951 | 3672 | 2951 | 722 |
| 1e-05 | 3083 | 2932 | 3346 | 2932 | 415 |
| 1e-04 | 3142 | 3039 | 3309 | 3039 | 270 |
| 1e-03 | 2997 | 2780 | 2904 | 2779 | 218 |
| 1e-02 | 2288 | 1951 | 1990 | 1951 | 337 |
| 1e-01 | 2653 | 2619 | 2630 | 2619 | 34 |
| 1e+00 | 3337 | 3331 | 3331 | 3331 | 6 |

### night: T [K]

| p [bar] | prodbin | sph | sphbeam | noang | max spread |
|---|---|---|---|---|---|
| 1e-06 | 1951 | 1427 | 1906 | 1434 | 524 |
| 1e-05 | 2109 | 1561 | 2066 | 1556 | 553 |
| 1e-04 | 2111 | 1629 | 1949 | 1624 | 486 |
| 1e-03 | 2028 | 1600 | 1616 | 1593 | 435 |
| 1e-02 | 2127 | 1768 | 1789 | 1767 | 360 |
| 1e-01 | 2661 | 2627 | 2631 | 2628 | 34 |
| 1e+00 | 3334 | 3329 | 3331 | 3329 | 5 |

### day - night contrast [K]

| p [bar] | prodbin | sph | sphbeam | noang |
|---|---|---|---|---|
| 1e-06 | 2370 | 2819 | 2373 | 2812 |
| 1e-05 | 2047 | 2498 | 2022 | 2503 |
| 1e-04 | 1716 | 2066 | 1757 | 2071 |
| 1e-03 | 1448 | 1797 | 1790 | 1804 |
| 1e-02 | 826 | 881 | 867 | 882 |
| 1e-01 | 102 | 107 | 106 | 106 |
| 1e+00 | 9 | 13 | 11 | 13 |

### differences from prodbin > 50 K

- sph, day, 1e-06 bar: -76 K
- noang, day, 1e-06 bar: -76 K
- sph, day, 1e-05 bar: -97 K
- sphbeam, day, 1e-05 bar: -69 K
- noang, day, 1e-05 bar: -97 K
- sph, day, 1e-04 bar: -131 K
- sphbeam, day, 1e-04 bar: -121 K
- noang, day, 1e-04 bar: -131 K
- sph, day, 1e-03 bar: -79 K
- sphbeam, day, 1e-03 bar: -70 K
- noang, day, 1e-03 bar: -79 K
- sph, day, 1e-02 bar: -305 K
- sphbeam, day, 1e-02 bar: -297 K
- noang, day, 1e-02 bar: -305 K
- sph, term, 1e-06 bar: -151 K
- sphbeam, term, 1e-06 bar: +570 K
- noang, term, 1e-06 bar: -152 K
- sph, term, 1e-05 bar: -151 K
- sphbeam, term, 1e-05 bar: +263 K
- noang, term, 1e-05 bar: -151 K
- sph, term, 1e-04 bar: -104 K
- sphbeam, term, 1e-04 bar: +167 K
- noang, term, 1e-04 bar: -104 K
- sph, term, 1e-03 bar: -217 K
- sphbeam, term, 1e-03 bar: -93 K
- noang, term, 1e-03 bar: -218 K
- sph, term, 1e-02 bar: -337 K
- sphbeam, term, 1e-02 bar: -298 K
- noang, term, 1e-02 bar: -337 K
- sph, night, 1e-06 bar: -524 K
- noang, night, 1e-06 bar: -518 K
- sph, night, 1e-05 bar: -548 K
- noang, night, 1e-05 bar: -553 K
- sph, night, 1e-04 bar: -481 K
- sphbeam, night, 1e-04 bar: -162 K
- noang, night, 1e-04 bar: -486 K
- sph, night, 1e-03 bar: -428 K
- sphbeam, night, 1e-03 bar: -412 K
- noang, night, 1e-03 bar: -435 K
- sph, night, 1e-02 bar: -359 K
- sphbeam, night, 1e-02 bar: -338 K
- noang, night, 1e-02 bar: -360 K

### horizontally averaged T [K] in the eta-capped shell (r/Rp 1.20-1.28)

| r/Rp | prodbin night | sph night | sphbeam night | noang night |
|---|---|---|---|---|
| 1.205 | 2713 | 2677 | 2682 | 2677 |
| 1.215 | 2528 | 2456 | 2466 | 2457 |
| 1.225 | 2327 | 2170 | 2182 | 2170 |
| 1.236 | 2167 | 1839 | 1857 | 1837 |
| 1.246 | 2105 | 1655 | 1685 | 1652 |
| 1.257 | 2045 | 1627 | 1641 | 1619 |
| 1.269 | 2063 | 1647 | 1742 | 1644 |
