# NOTE 2026-11-16 Orion -> viper (USER): the clean scheme moves back to viper, in three phases

USER: "slowly move the clean scheme writing and testing to viper, since the basic scheme is there and gpu testing
is important now". Full handover: docs/handover/HANDOVER-2026-11-16-clean-scheme-orion-to-viper.md (scheme, key
names and defaults, binding rules, ranked open problems, test harnesses, paths).

1. Phase 1 (now): viper runs the GPU tests (TASK cc_gpu1 / -a100, stage 1) and the sp/cs GPU divergence
   debugging (TASK cc_gpu_debug 8e010e76d, D1-D4 on 094488c1f; origin/clean-cc is now push 12 0eeac5fdb = + the WB limiter fix, sp/cs GPU paths unchanged). Orion still owns clean-cc and
   audits the sp/cs device paths in parallel (clean-cc-gpufix). The decisions TASK (cc_gpu_decisions) stays for
   Raven later.
2. Phase 2: Orion finishes its in-flight items (time-order diagnosis, GPU fix if the audit
   finds it, clump test merge; time-dependent rays awaits the user), pushes them as origin/clean-cc after a green
   regression, and packages the regression / restart / rw / IDORT-page harnesses on bsg-files-1009. viper reads
   the code and reruns the regression once on Raven CPU.
3. Phase 3: after Orion's "phase 3" NOTE (last hash, unmerged branches, state of each open problem) viper owns
   clean-cc development and testing; Orion supports on CPU on request.

Protocol unchanged: NOTE/TASK files on bsg-files-1009; one published branch origin/clean-cc, advanced only after
a green regression and only by the current owner (Orion until phase 3); restyle through origin/clean-1114, merged
into clean-cc by the clean owner after each batch (fd5106e4b). origin/clean-cc-int (4680ce500) stays frozen.
