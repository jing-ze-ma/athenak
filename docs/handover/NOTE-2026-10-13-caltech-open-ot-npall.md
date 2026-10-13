# NOTE 2026-10-13 from Caltech: open-ot NPALL (re TASK-2026-10-13-caltech-open-ot-npall.md)

## QUEUED (10-13 ~00:55 PDT)
- Binary: 281dd686 + patch/0001-0004 (git apply; no kokkos change in these) -> athena_cpu_omp_281dd686_npall4; md5 follows when the build ends.
- CVALL FULL step stopped as asked (runner returns after BIT); the remaining 18 BIT dirs (14 re-dos + 4) run in CVALL round 4429437; extra CVALL rounds cancelled.
- NPALL: debug-QoS chain 4434594 -> 4434595 -> 4434596 -> 4434597 -> 4434598 -> 4434599 (after CVALL round 2), P = 10, 32 cores, finished dirs skipped.
  Order = your priority: sp/*_T, lad/*, other sp, cart. Skipped as long: cart/cloud/FIN3_r12g2, cart/ddf/FIN3_d, cart/ddiff/FIN3_d (pipe is not in list.txt).
- Caltech change to run_npall.sh: the shock-ref sed also maps /viper/ptmp2/jinma/fin3_1012/shockref/ (only the Raven path was mapped; 3 of the rs decks use the Viper path).
