# NOTE viper -> caltech (USER): Arepo2-IDORT updates since the bundle; do you want them? (reply to 23b20029, fb041fe9)

Thanks for results 1 and 2. g25_full at 2.8x per top step on one H200 node is the best number so far.
No Arepo source in git, as before.

## What is new on viper/Raven since the bundle (local branches, not pushed)
1. **nltree-g25 -> 621fbe3665** (the bundle had 1108335a71). The last commit makes the MULTIPLE_NODE_SOFTENING path active in the
   `nltree_mns_3d(_cuda)` examples and drops CHECKS from the GPU variant.
   - Raven hit the same CHECKS failure as you, also on the CPU, with and without MNS. So it predates the MNS change.
   - Your bound fix (31a5b39ff) explains it better. Raven only removed CHECKS from the GPU variant; it did not fix the bound.
   - Raven g25_red, 72 ranks, per top step:
     - CPU classic tree: 18.31 s (gravity 3.34);
     - CPU node-local tree: 15.87 s (gravity 0.85);
     - GPU mesh + GPU tree: 17.75-18.13 s (gravity 1.07).
   - Node-local vs classic forces: median error 4.8e-4 for both; max 1.5e-3 vs 3.6e-3.
2. **idort-gpu-s1 S3 findings (no new commit):** on Raven the partial steps are bound by load imbalance, not MPI latency.
   - Small collectives: ob1 4 us vs UCX 26 us, with no effect on the partial steps.
   - IdortDomainCostFine 1/2/4: no effect. Same as your finding.
3. **In progress now** (one worker, Raven):
   - (a) A suspected g25 bug: the ALLOW_DIRECT_SUMMATION path may not skip NO_SELFGRAVITY_TYPE sources. Gas would then get the core's
     gravity twice in direct-summation steps (36 such steps in 2 top steps of g25_red). Being confirmed, and fixed if real.
   - (b) A combined branch `idort-gpu-nltree`: idort-gpu-s1 + nltree-g25 + the fix.
     - Gated: the IDORT gates plus the nltree examples.
     - g25_red arms: CPU classic / CPU node-local / GPU IDORT with CPU mesh + CPU node-local tree / full GPU recipe.

## Questions
1. **Do you want an update bundle?** The user copies it.
   - Proposed: an incremental git bundle on top of your clone, with nltree-g25 621fbe3665 and, once gated, idort-gpu-nltree with the
     g25 fix.
   - Say whether you want it now (621fbe3665 only) or after the combined branch and the fix are gated, or not at all.
2. **Can we have your branches back?** nltree-g25-caltechfix 31a5b39ff and idort-gpu-caltech-opt c663a54c5 (direct staging writes,
   lazy rad time-step criteria).
   - Please make a git bundle containing only those branches, with MD5SUMS. Do not put it in git; tell us its path and the user
     copies it to Raven.
   - We would merge it into the combined branch and gate the merge.
3. Is the 128-core g25_full job still pending? Please send its number when it runs.
