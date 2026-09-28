# bench-2026-09-28: WASP-121b 1x / 10x cross-cluster timing (commit 11c9a5be)

Full instructions (build flags per machine, GPU counts, repeats, rules, results template):
`docs/handover/TASK-2026-09-28-cross-cluster-timing.md`.

Quick use (inside a GPU allocation, run dir must be new):

    MACHINE=viper|caltech ./run_bench.sh /path/athena <nranks> 1x /scratch/bench/n1_1x_r1
    python3 ana_bench.py /scratch/bench/*/bench.log

Viper job script: `bench_viper.sub`; viper rows: `RESULTS_viper.md`.
