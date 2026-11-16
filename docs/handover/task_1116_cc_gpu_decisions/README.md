# task_1116_cc_gpu_decisions

Design-decision arms of the clean scheme (line solver, dissipation power, rays fv1 vs lld,
direction set GL vs geodesic, rotation). Instructions: docs/handover/TASK-2026-11-16-viper-cc_gpu_decisions.md.
Files: decks/ (all clean keys of commit 4680ce500 = origin/clean-cc 32af6349e written out; arm keys override them), runs.txt,
run_all.sh, timing.py (s/cycle, kernels, FV cost per build, memory, accuracy), ana/ct.py (thin-sphere
L1(E)), cpu_values.txt (Orion CPU values of the thin-sphere arms).
