# NOTE 2026-11-21gy orion -> caltech: make CUDA 12.9 builds fast (USER: "make build fast")

1. Give CUDA 12.9 builds their OWN ccache: CCACHE_DIR=/resnick/groups/carnegie_poc/jingze/ccache_cu129 (new, empty), keep
   CCACHE_BASEDIR/NOHASHDIR as before. Never point a 12.9 build at the old shared cache (12.2 objects). The compiler path is in
   the ccache hash, but a separate dir removes any doubt and keeps the hit rate high for 12.9 trees.
2. Build on a whole free node: srun -N1 -c 32 (the holders have 32 CPUs) + make -j32; one target per node, several in parallel.
3. Right now, in parallel with the BV runs (NOTE gx), pre-build with CUDA 12.9 into separate dirs so the next launch doesn't wait:
   - B_v204n: v20.4 code head 7cbc13437 (he_star_m1, same flags as B_v203n) — for later, not for production yet.
   Report md5 + build seconds; then run the tiny BSG smoke gate on it (pass 1 must be 4.020e-06).
Every future GPU binary must pass that 1-minute smoke gate before use.

-- Orion coordinator
