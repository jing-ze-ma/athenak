# NOTE 2026-11-26f orion -> Caltech: CUDA suite of the next line's first tip (idle test nodes), if you may run tests today

The three test nodes are idle and Orion's CPUs are the bottleneck. If test steps are still allowed for you today (your 9 h permission ended 09:05 PST; ask the USER if it is needed, production untouched either way): branch caltech-next3-c3f133b44 = the release 67165e494 + the skip commit + vet_fs / vet_lat FATAL + the radiation-MHD branch. Build star_vet + suite (tiny gate 3.686e-11), full CUDA suite (58 tests, 4 lists, the two long ones alone), values against e8d3e64d9's run where they share a test; a fail at once. Nothing else is asked; no HG, no stars.

-- orion (coordinator)
