## Part 1: optically thin RE roots (s = stable, u = unstable), T grid 300-4000 K

Fixed P (isobaric, Field). Columns: P = 1e-6, 1e-4, 1e-2 dyn/cm^2 (rho ~ 1e-17..1e-13 at these T); fixed-rho roots are listed in part1_roots_rho.txt.

| star | mode | r/R | W | P=1e-6 | P=1e-4 | P=1e-2 | two stable (Tw/Tc) at 1e-4 | rho_c/rho_w |
|---|---|---|---|---|---|---|---|---|
| golden16 | clamp | 1.5 | 0.1273 | 3358s | 3358s | 3358s | no | - |
| golden16 | clamp | 2.0 | 0.0670 | 2940s | 2940s | 2940s | no | - |
| golden16 | clamp | 3.0 | 0.0286 | 939s,1089u,1122s,1339u,2744s | 939s,1089u,1122s,1339u,2744s | 939s,1089u,1122s,1339u,2744s | 2744 / 939 | 5.4 |
| golden16 | clamp | 4.0 | 0.0159 | 795s,1394u,1494s,1519u,2579s | 795s,1394u,1494s,1519u,2579s | 795s,1394u,1494s,1519u,2579s | 2579 / 795 | 6.0 |
| golden16 | clamp | 5.0 | 0.0101 | 708s,1757u,1928s | 708s,1757u,1928s | 708s,1757u,1928s | 1928 / 708 | 5.0 |
| golden16 | extrap | 1.5 | 0.1273 | 2672s,3060u,3377s | 2783s,3025u,3377s | 3358s | 3377 / 2783 | 1.2 |
| golden16 | extrap | 2.0 | 0.0670 | 2617s,3133u,3184s | 2708s | 2940s | no | - |
| golden16 | extrap | 3.0 | 0.0286 | 941s,1089u,1122s,1344u,2552s | 940s,1089u,1122s,1342u,2596s | 939s,1089u,1122s,1339u,2744s | 2596 / 940 | 5.1 |
| golden16 | extrap | 4.0 | 0.0159 | 795s,1396u,1422s,1503u,2509s | 795s,1395u,1461s,1504u,2521s | 795s,1394u,1494s,1519u,2579s | 2521 / 795 | 5.9 |
| golden16 | extrap | 5.0 | 0.0101 | 708s,1548u,1879s | 708s,1591u,1891s | 708s,1757u,1928s | 1891 / 708 | 4.9 |
| m20lgl5.5 | clamp | 1.5 | 0.1273 | 2890s | 2890s | 2890s | no | - |
| m20lgl5.5 | clamp | 2.0 | 0.0670 | 2753s | 2753s | 2753s | no | - |
| m20lgl5.5 | clamp | 3.0 | 0.0286 | 894s,1385u,2038s,2191u,2496s | 894s,1385u,2038s,2191u,2496s | 894s,1385u,2038s,2191u,2496s | 2496 / 894 | 5.2 |
| m20lgl5.5 | clamp | 4.0 | 0.0159 | 763s,1797u,1902s | 763s,1797u,1902s | 763s,1797u,1902s | 1902 / 763 | 4.6 |
| m20lgl5.5 | clamp | 5.0 | 0.0101 | 684s | 684s | 684s | no | - |
| m20lgl5.5 | extrap | 1.5 | 0.1273 | 2610s | 2696s | 2890s | no | - |
| m20lgl5.5 | extrap | 2.0 | 0.0670 | 2557s | 2603s | 2753s | no | - |
| m20lgl5.5 | extrap | 3.0 | 0.0286 | 895s,1386u,1911s,2157u,2453s | 895s,1385u,1940s,2171u,2464s | 894s,1385u,2038s,2191u,2496s | 2464 / 895 | 5.1 |
| m20lgl5.5 | extrap | 4.0 | 0.0159 | 763s,1555u,1856s | 763s,1606u,1856s | 763s,1797u,1902s | 1856 / 763 | 4.4 |
| m20lgl5.5 | extrap | 5.0 | 0.0101 | 684s,1603u,1803s | 684s,1707u,1753s | 684s | 1753 / 684 | 4.3 |

## Part 2: gas only, max over (case, mode, grad, f, T_c choice)

| model | r/R | Tbar | rho_bar | Gamma_w(Sob) max | Gamma_c gas(Sob) max | Gamma_bar gas max | Gamma_c thin max | RE pair (mft clamp, f=0.1) |
|---|---|---|---|---|---|---|---|---|
| m15lgl5.1 | 1.5 | 2449 | 6.4e-13 | 0.004 | 0.045 | 0.028 | 0.218 | none |
| m15lgl5.1 | 2.0 | 2194 | 8.1e-14 | 0.011 | 0.089 | 0.053 | 0.216 | none |
| m15lgl5.1 | 3.0 | 2083 | 1.0e-14 | 0.032 | 0.136 | 0.087 | 0.216 | 2744/939 K, f_m 0.38 |
| m15lgl5.1 | 4.0 | 1552 | 3.3e-15 | 0.196 | 0.169 | 0.192 | 0.216 | 2579/795 K, f_m 0.40 |
| m15lgl5.1 | 5.0 | 1086 | 1.5e-15 | 0.158 | 0.016 | 0.130 | 0.018 | 1928/708 K, f_m 0.36 |
| m20lgl5.3 | 1.5 | 2438 | 5.3e-13 | 0.005 | 0.050 | 0.030 | 0.251 | none |
| m20lgl5.3 | 2.0 | 2185 | 6.6e-14 | 0.012 | 0.099 | 0.059 | 0.248 | none |
| m20lgl5.3 | 3.0 | 2078 | 8.3e-15 | 0.036 | 0.153 | 0.098 | 0.248 | 2735/936 K, f_m 0.38 |
| m20lgl5.3 | 4.0 | 1519 | 2.6e-15 | 0.220 | 0.193 | 0.216 | 0.248 | 2564/793 K, f_m 0.40 |
| m20lgl5.3 | 5.0 | 1079 | 1.2e-15 | 0.187 | 0.019 | 0.156 | 0.021 | 1921/706 K, f_m 0.36 |
| m20lgl5.5 | 1.5 | 2339 | 1.6e-12 | 0.007 | 0.072 | 0.043 | 0.414 | none |
| m20lgl5.5 | 2.0 | 2094 | 2.7e-13 | 0.034 | 0.149 | 0.088 | 0.411 | none |
| m20lgl5.5 | 3.0 | 1549 | 4.5e-14 | 0.144 | 0.241 | 0.168 | 0.409 | 2496/894 K, f_m 0.36 |
| m20lgl5.5 | 4.0 | 1280 | 1.6e-14 | 0.307 | 0.033 | 0.274 | 0.107 | 1902/763 K, f_m 0.34 |
| m20lgl5.5 | 5.0 | 1123 | 7.9e-15 | 0.260 | 0.027 | 0.248 | 0.044 | none |

## Part 3a: where the condensation criteria are met (T_cond optimistic / nominal)

| model | r/R | T_d p=1 | T_d p=0 | T_d p=-1 | Al2O3 (1600/1400) grain p=1, 0, -1 | silicate (1200/1000) grain p=1, 0, -1 |
|---|---|---|---|---|---|---|
| m15lgl5.1 | 1.5 | 2719 | 2453 | 2066 | no, no, no | no, no, no |
| m15lgl5.1 | 2.0 | 2391 | 2089 | 1668 | no, no, no | no, no, no |
| m15lgl5.1 | 3.0 | 2017 | 1689 | 1256 | no, no, yes | no, no, no |
| m15lgl5.1 | 4.0 | 1793 | 1458 | 1032 | no, opt, yes | no, no, opt |
| m15lgl5.1 | 5.0 | 1638 | 1302 | 888 | no, yes, yes | no, no, yes |
| m20lgl5.3 | 1.5 | 2707 | 2442 | 2056 | no, no, no | no, no, no |
| m20lgl5.3 | 2.0 | 2380 | 2080 | 1660 | no, no, no | no, no, no |
| m20lgl5.3 | 3.0 | 2008 | 1681 | 1250 | no, no, yes | no, no, no |
| m20lgl5.3 | 4.0 | 1785 | 1451 | 1027 | no, opt, yes | no, no, opt |
| m20lgl5.3 | 5.0 | 1631 | 1296 | 884 | no, yes, yes | no, no, yes |
| m20lgl5.5 | 1.5 | 2501 | 2256 | 1900 | no, no, no | no, no, no |
| m20lgl5.5 | 2.0 | 2199 | 1921 | 1534 | no, no, opt | no, no, no |
| m20lgl5.5 | 3.0 | 1855 | 1553 | 1155 | no, opt, yes | no, no, opt |
| m20lgl5.5 | 4.0 | 1649 | 1341 | 949 | no, yes, yes | no, no, yes |
| m20lgl5.5 | 5.0 | 1506 | 1197 | 816 | opt, yes, yes | no, opt, yes |

(yes = below the nominal (lower) T_cond; opt = only below the optimistic (upper) T_cond.)  Gas criterion: T_c (cool-phase gas T) < T_cond; with the RE pair T_c = 700-940 K at 3-5 R, so both species pass the gas criterion wherever a two-phase RE solution exists; in the T_c scan, 1000/1200/1400 K.


## Part 3b: max Gamma_c and Gamma_bar (T_cond optimistic; porosity on, l = 0.01 r or 0.1 r; l=0 column = no porosity)

| model | species | criterion | r range | max Gamma_c | where | max Gamma_bar | where | max Gamma_c, l=0 | max Gamma_bar, l=0 |
|---|---|---|---|---|---|---|---|---|---|
| m15lgl5.1 | Al2O3 | gas | 2-4 | 0.342 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.243 | r=4.0 mft12/extrap v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.342 | 0.243 |
| m15lgl5.1 | Al2O3 | gas | 1.5-5 | 0.342 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.243 | r=4.0 mft12/extrap v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.342 | 0.243 |
| m15lgl5.1 | Al2O3 | grain_p1 | 2-4 | 0.169 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.192 | r=4.0 mft12/clamp v/H f=0.01 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.169 | 0.192 |
| m15lgl5.1 | Al2O3 | grain_p1 | 1.5-5 | 0.169 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.192 | r=4.0 mft12/clamp v/H f=0.01 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.169 | 0.192 |
| m15lgl5.1 | Al2O3 | grain_p0 | 2-4 | 0.342 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.243 | r=4.0 mft12/extrap v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.342 | 0.243 |
| m15lgl5.1 | Al2O3 | grain_p0 | 1.5-5 | 0.342 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.243 | r=4.0 mft12/extrap v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.342 | 0.243 |
| m15lgl5.1 | Al2O3 | grain_pm1 | 2-4 | 0.342 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.243 | r=4.0 mft12/extrap v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.342 | 0.243 |
| m15lgl5.1 | Al2O3 | grain_pm1 | 1.5-5 | 0.342 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.243 | r=4.0 mft12/extrap v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.342 | 0.243 |
| m15lgl5.1 | silicate | gas | 2-4 | 6.797 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1000 fc=1.0 kd=3000 l=0.01 | 4.830 | r=4.0 mft12/clamp v/H f=0.3 RE Tc=795 fc=1.0 kd=3000 l=0.01 | 6.933 | 4.998 |
| m15lgl5.1 | silicate | gas | 1.5-5 | 6.880 | r=5.0 mft12/extrap v/H f=0.3 scan Tc=1000 fc=1.0 kd=3000 l=0.01 | 4.830 | r=4.0 mft12/clamp v/H f=0.3 RE Tc=795 fc=1.0 kd=3000 l=0.01 | 6.935 | 4.998 |
| m15lgl5.1 | silicate | grain_p1 | 2-4 | 0.169 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.192 | r=4.0 mft12/clamp v/H f=0.01 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.169 | 0.192 |
| m15lgl5.1 | silicate | grain_p1 | 1.5-5 | 0.169 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.192 | r=4.0 mft12/clamp v/H f=0.01 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.169 | 0.192 |
| m15lgl5.1 | silicate | grain_p0 | 2-4 | 0.169 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.192 | r=4.0 mft12/clamp v/H f=0.01 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.169 | 0.192 |
| m15lgl5.1 | silicate | grain_p0 | 1.5-5 | 0.169 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.192 | r=4.0 mft12/clamp v/H f=0.01 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.169 | 0.192 |
| m15lgl5.1 | silicate | grain_pm1 | 2-4 | 6.974 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 4.830 | r=4.0 mft12/clamp v/H f=0.3 RE Tc=795 fc=1.0 kd=3000 l=0.01 | 7.088 | 4.998 |
| m15lgl5.1 | silicate | grain_pm1 | 1.5-5 | 6.974 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 4.830 | r=4.0 mft12/clamp v/H f=0.3 RE Tc=795 fc=1.0 kd=3000 l=0.01 | 7.088 | 4.998 |
| m20lgl5.3 | Al2O3 | gas | 2-4 | 0.392 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.274 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.393 | 0.274 |
| m20lgl5.3 | Al2O3 | gas | 1.5-5 | 0.392 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.274 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.393 | 0.274 |
| m20lgl5.3 | Al2O3 | grain_p1 | 2-4 | 0.193 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.216 | r=4.0 mft12/clamp v/H f=0.01 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.193 | 0.216 |
| m20lgl5.3 | Al2O3 | grain_p1 | 1.5-5 | 0.193 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.216 | r=4.0 mft12/clamp v/H f=0.01 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.193 | 0.216 |
| m20lgl5.3 | Al2O3 | grain_p0 | 2-4 | 0.392 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.274 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.393 | 0.274 |
| m20lgl5.3 | Al2O3 | grain_p0 | 1.5-5 | 0.392 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.274 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.393 | 0.274 |
| m20lgl5.3 | Al2O3 | grain_pm1 | 2-4 | 0.392 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.274 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.393 | 0.274 |
| m20lgl5.3 | Al2O3 | grain_pm1 | 1.5-5 | 0.392 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.274 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.393 | 0.274 |
| m20lgl5.3 | silicate | gas | 2-4 | 7.848 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1000 fc=1.0 kd=3000 l=0.01 | 5.557 | r=4.0 mft12/clamp v/H f=0.3 RE Tc=793 fc=1.0 kd=3000 l=0.01 | 8.024 | 5.779 |
| m20lgl5.3 | silicate | gas | 1.5-5 | 7.954 | r=5.0 mft12/extrap v/H f=0.3 scan Tc=1000 fc=1.0 kd=3000 l=0.01 | 5.557 | r=4.0 mft12/clamp v/H f=0.3 RE Tc=793 fc=1.0 kd=3000 l=0.01 | 8.027 | 5.779 |
| m20lgl5.3 | silicate | grain_p1 | 2-4 | 0.193 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.216 | r=4.0 mft12/clamp v/H f=0.01 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.193 | 0.216 |
| m20lgl5.3 | silicate | grain_p1 | 1.5-5 | 0.193 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.216 | r=4.0 mft12/clamp v/H f=0.01 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.193 | 0.216 |
| m20lgl5.3 | silicate | grain_p0 | 2-4 | 0.193 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.216 | r=4.0 mft12/clamp v/H f=0.01 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.193 | 0.216 |
| m20lgl5.3 | silicate | grain_p0 | 1.5-5 | 0.193 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.216 | r=4.0 mft12/clamp v/H f=0.01 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.193 | 0.216 |
| m20lgl5.3 | silicate | grain_pm1 | 2-4 | 8.054 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 5.557 | r=4.0 mft12/clamp v/H f=0.3 RE Tc=793 fc=1.0 kd=3000 l=0.01 | 8.201 | 5.779 |
| m20lgl5.3 | silicate | grain_pm1 | 1.5-5 | 8.054 | r=4.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 5.557 | r=4.0 mft12/clamp v/H f=0.3 RE Tc=793 fc=1.0 kd=3000 l=0.01 | 8.201 | 5.779 |
| m20lgl5.5 | Al2O3 | gas | 2-4 | 0.612 | r=3.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.365 | r=4.0 mft12/clamp v/H f=0.3 RE Tc=763 fc=1.0 kd=3000 l=0.01 | 0.614 | 0.365 |
| m20lgl5.5 | Al2O3 | gas | 1.5-5 | 0.612 | r=3.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.365 | r=4.0 mft12/clamp v/H f=0.3 RE Tc=763 fc=1.0 kd=3000 l=0.01 | 0.614 | 0.365 |
| m20lgl5.5 | Al2O3 | grain_p1 | 2-4 | 0.240 | r=3.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.274 | r=4.0 mft12/clamp v/H f=0.01 RE Tc=763 fc=0.1 kd=300 l=0.01 | 0.241 | 0.274 |
| m20lgl5.5 | Al2O3 | grain_p1 | 1.5-5 | 0.400 | r=5.0 mft12/extrap v/H f=0.3 RE Tc=684 fc=1.0 kd=3000 l=0.01 | 0.352 | r=5.0 mft12/extrap v/H f=0.3 RE Tc=684 fc=1.0 kd=3000 l=0.01 | 0.400 | 0.352 |
| m20lgl5.5 | Al2O3 | grain_p0 | 2-4 | 0.612 | r=3.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.365 | r=4.0 mft12/clamp v/H f=0.3 RE Tc=763 fc=1.0 kd=3000 l=0.01 | 0.614 | 0.365 |
| m20lgl5.5 | Al2O3 | grain_p0 | 1.5-5 | 0.612 | r=3.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.365 | r=4.0 mft12/clamp v/H f=0.3 RE Tc=763 fc=1.0 kd=3000 l=0.01 | 0.614 | 0.365 |
| m20lgl5.5 | Al2O3 | grain_pm1 | 2-4 | 0.612 | r=3.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.365 | r=4.0 mft12/clamp v/H f=0.3 RE Tc=763 fc=1.0 kd=3000 l=0.01 | 0.614 | 0.365 |
| m20lgl5.5 | Al2O3 | grain_pm1 | 1.5-5 | 0.612 | r=3.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=1.0 kd=3000 l=0.01 | 0.365 | r=4.0 mft12/clamp v/H f=0.3 RE Tc=763 fc=1.0 kd=3000 l=0.01 | 0.614 | 0.365 |
| m20lgl5.5 | silicate | gas | 2-4 | 14.692 | r=4.0 mft12/extrap v/H f=0.3 scan Tc=1000 fc=1.0 kd=3000 l=0.01 | 9.623 | r=4.0 mft12/clamp v/H f=0.3 RE Tc=763 fc=1.0 kd=3000 l=0.01 | 14.979 | 10.316 |
| m20lgl5.5 | silicate | gas | 1.5-5 | 14.692 | r=4.0 mft12/extrap v/H f=0.3 scan Tc=1000 fc=1.0 kd=3000 l=0.01 | 9.726 | r=5.0 mft12/extrap v/H f=0.3 RE Tc=684 fc=1.0 kd=3000 l=0.01 | 14.979 | 10.316 |
| m20lgl5.5 | silicate | grain_p1 | 2-4 | 0.240 | r=3.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.274 | r=4.0 mft12/clamp v/H f=0.01 RE Tc=763 fc=0.1 kd=300 l=0.01 | 0.241 | 0.274 |
| m20lgl5.5 | silicate | grain_p1 | 1.5-5 | 0.240 | r=3.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.274 | r=4.0 mft12/clamp v/H f=0.01 RE Tc=763 fc=0.1 kd=300 l=0.01 | 0.241 | 0.274 |
| m20lgl5.5 | silicate | grain_p0 | 2-4 | 0.240 | r=3.0 mft12/clamp v/H f=0.3 scan Tc=1400 fc=0.1 kd=300 l=0.01 | 0.274 | r=4.0 mft12/clamp v/H f=0.01 RE Tc=763 fc=0.1 kd=300 l=0.01 | 0.241 | 0.274 |
| m20lgl5.5 | silicate | grain_p0 | 1.5-5 | 14.656 | r=5.0 mft12/extrap v/H f=0.3 RE Tc=684 fc=1.0 kd=3000 l=0.01 | 9.726 | r=5.0 mft12/extrap v/H f=0.3 RE Tc=684 fc=1.0 kd=3000 l=0.01 | 14.973 | 9.934 |
| m20lgl5.5 | silicate | grain_pm1 | 2-4 | 14.692 | r=4.0 mft12/extrap v/H f=0.3 scan Tc=1000 fc=1.0 kd=3000 l=0.01 | 9.623 | r=4.0 mft12/clamp v/H f=0.3 RE Tc=763 fc=1.0 kd=3000 l=0.01 | 15.187 | 10.316 |
| m20lgl5.5 | silicate | grain_pm1 | 1.5-5 | 14.692 | r=4.0 mft12/extrap v/H f=0.3 scan Tc=1000 fc=1.0 kd=3000 l=0.01 | 9.726 | r=5.0 mft12/extrap v/H f=0.3 RE Tc=684 fc=1.0 kd=3000 l=0.01 | 15.187 | 10.316 |

## Part 3c: minimum delta*f_cond*kappa_d [cm^2/g gas] for Gamma_c > 1 (case mft, clamp, grad v/H, f = 0.1; cool phase = RE pair if it exists, else T_c = 1200 K scan); inf = impossible: max effective opacity 1/(rho_c l) < kappa_Edd

| model | r/R | kappa_Edd | T_c | rho_c | kappa_gas,c(Sob) | X_req l=0 | X_req l=0.01r | 1/(rho_c l) l=0.01r | X_req l=0.1r | 1/(rho_c l) l=0.1r |
|---|---|---|---|---|---|---|---|---|---|---|
| m15lgl5.1 | 1.5 | 1.73 | 1200 | 1.9e-12 | 0.009 | 1.73 | inf | 1 | inf | 0.1 |
| m15lgl5.1 | 2.0 | 1.73 | 1200 | 2.2e-13 | 0.016 | 1.72 | 2.14 | 5 | inf | 0.5 |
| m15lgl5.1 | 3.0 | 1.73 | 939 | 3.9e-14 | 0.010 | 1.72 | 1.81 | 18 | 5.37 | 1.8 |
| m15lgl5.1 | 4.0 | 1.73 | 795 | 1.3e-14 | 0.011 | 1.72 | 1.76 | 40 | 2.26 | 4.0 |
| m15lgl5.1 | 5.0 | 1.73 | 708 | 5.4e-15 | 0.013 | 1.72 | 1.74 | 79 | 1.95 | 7.9 |
| m20lgl5.3 | 1.5 | 1.50 | 1200 | 1.6e-12 | 0.009 | 1.49 | inf | 1 | inf | 0.1 |
| m20lgl5.3 | 2.0 | 1.50 | 1200 | 1.8e-13 | 0.016 | 1.48 | 1.79 | 5 | inf | 0.5 |
| m20lgl5.3 | 3.0 | 1.50 | 936 | 3.1e-14 | 0.011 | 1.49 | 1.55 | 18 | 3.13 | 1.8 |
| m20lgl5.3 | 4.0 | 1.50 | 793 | 1.1e-14 | 0.011 | 1.49 | 1.52 | 41 | 1.86 | 4.1 |
| m20lgl5.3 | 5.0 | 1.50 | 707 | 4.3e-15 | 0.014 | 1.48 | 1.50 | 80 | 1.65 | 8.0 |
| m20lgl5.5 | 1.5 | 0.80 | 1200 | 4.5e-12 | 0.004 | 0.80 | inf | 0 | inf | 0.0 |
| m20lgl5.5 | 2.0 | 0.80 | 1200 | 7.2e-13 | 0.008 | 0.80 | inf | 1 | inf | 0.1 |
| m20lgl5.5 | 3.0 | 0.80 | 895 | 1.6e-13 | 0.005 | 0.80 | 0.99 | 2 | inf | 0.2 |
| m20lgl5.5 | 4.0 | 0.80 | 763 | 5.4e-14 | 0.005 | 0.80 | 0.87 | 5 | inf | 0.5 |

Supply (delta_max * f_cond=1 * kappa_d=3000): Al2O3 1e-4*3000 = 0.30; silicate 4e-3*3000 = 12 cm^2/g gas (kappa_d = 300: 0.03 / 1.2).

