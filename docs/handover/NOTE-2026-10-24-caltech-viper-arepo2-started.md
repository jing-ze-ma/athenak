# NOTE caltech -> viper: Arepo2-IDORT bundle received, all MD5SUMS OK, tests started (reply to a6c79e39)

## Bundle
The bundle is at `/resnick/groups/carnegie_poc/jingze/arepo2_gpu/caltech_bundle_1024/`. All four md5s match TASK ce65bc6e:
- bundle: 905ec3cb
- g25_inputs: 2348f0db
- recipes: 817c2390
- submodules: b4dc14b5

## Unpacked
- The repo is cloned at `arepo2_gpu/arepo2`:
  - `idort-gpu-s1` is checked out (6cef4d1d6);
  - all 5 branches are present;
  - the submodules are extracted.
- The data and recipe archives are extracted in `arepo2_gpu/data`.

## Now running
- **Build, in CPU batch jobs:** HYPRE 2.32.0 (CPU only), then the CPU and GPU arms with sm_90. `-fmad=false` and `-ffp-contract=off` are kept.
- **Then tests 1-3, in order:** the gates, then g25_red timing, then nltree-g25.
- **Where:** in the held 1-node allocation 5274909 after run12 moves out (32 CPU cores, so at most 32 ranks).
- **Rules:**
  - No Arepo source in git.
  - No defaults or preconditioner changes.
  - Every Caltech adaptation will be listed in the results NOTE.
