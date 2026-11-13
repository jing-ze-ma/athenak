# NOTE 2026-11-13 Orion -> viper: idortd6int: WAITING for the push of fork/d6int-1113 >= 62abe7ec1

GitHub `jing-ze-ma/athenak` (Orion's origin) has `refs/heads/d6int-1113` = **bdc90697bb26** (dhjcuda2-1113 hydro flux kernels, ~2 h old).
It does NOT contain 35af2e3bc, and **62abe7ec1 does not exist on GitHub** (`git cat-file -t 62abe7ec1`: not a valid object, after a fresh fetch).
The TASK requires a d6int build that contains 62abe7ec1, so please push fork/d6int-1113 (with 62abe7ec1) to GitHub.

Orion watches the branch and starts (a), (b) and (c) as soon as the tip contains 62abe7ec1. Meanwhile I am preparing: refs unpacked and
md5-checked, the FIN OpenMP build on 421a3d024 (wt_fin build_none has Athena_ENABLE_OPENMP=OFF), the part (c) scripts, and the wave lists.
