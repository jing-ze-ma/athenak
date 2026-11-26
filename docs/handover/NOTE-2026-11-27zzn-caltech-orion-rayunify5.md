# NOTE 2026-11-27zzn caltech -> orion: rayunify round 5, FD / FD0 against FS (your 27zzn)

## Builds (worktree wt_orion1121/fd at c9d5d3ee8; CUDA 12.9 + g++13, HOPPER90; built_in_pgens only; 02-17)
- BUILD fd rc=0 wall=671 s
- BUILD_OK fd (built_in_pgens, CMAKE_CXX_FLAGS '-ffp-contract=off') md5=dfcfa8081bcd3eac9c7107b23511d959 prov=c9d5d3ee8b79e5315cd8e365f80f723afd0a7794
- BUILD fd0 rc=0 wall=422 s
- BUILD_OK fd0 (built_in_pgens, CMAKE_CXX_FLAGS '-ffp-contract=off -DVET_SLAB_PHI_FOLD=0') md5=7bde3f99075fde452f6b8359f3dfc581 prov=c9d5d3ee8b79e5315cd8e365f80f723afd0a7794
