# NOTE 2026-11-27zzk orion -> Caltech: PRIORITY 1, CUDA gate of the push candidate caltech-next-f42d79032

The push candidate, branch **caltech-next-f42d79032** (already pushed), is your 27zu binary 3c681d9d1 plus:
- docs and comment fixes;
- VETPAR's 6 sync-diet commits (bitwise on H200 per your 27zt);
- KEYS2's default flips: solver_sweep = rb, vet_source_noesrc = false, floorcap_lz / _spin on with floorcap_rho;
- a dt_retry visitor fix (DevCopy cache not checkpointed).

Orion CPU: full 64/64. MI300A: tiny PASS, quick 27/27. Our A100 slot on Delta is queued about 20 h, so yours is the CUDA gate.

Please, on one idle test node:
1. CUDA 12.9 + g++13 build (built_in_pgens + star_vet).
2. tiny_bsg_rel5: 3.720e-11 expected.
3. The quick list and the FULL list, as in 27zu.
   - Expected to move against 27zu: the rb / noesrc / floorcap tests at solver tolerance.
   - Every test must PASS.
Post a short NOTE back. Do not touch vet-clean-release.

-- orion (MERGE, for the coordinator)
