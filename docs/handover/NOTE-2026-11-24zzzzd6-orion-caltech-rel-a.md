# NOTE 2026-11-24zzzzd6 orion -> Caltech: GPU gate of the merged code fixes (replaces zzzzb5's suite)

Branch `caltech-rel-a-c880e2703` = GitHub release a5a4a67d3 + sdirk2 fix (7 commits) + cs thin-flux fix + ray-hole cure + the beams pin + docs (Orion CPU: full suite 38/38; key-off 258 files bitwise to 95813cf85 and to the release). **It replaces the sdorder5 GPU suite of zzzzb5** (skip that if it hasn't started; if running, let it end; a build done there is fine to keep).
1. CUDA 12.9 + g++ 13, own build dir: built_in suite binary, star_vet binary; he_star_vet tiny gate 3.686e-11; md5 + .prov.
2. Full GPU suite (full list, 38 tests) on one idle node: n/N, names of fails, every value that moved against the pin, CLEAN WARNING count. The idort_beams pin is the new one.
3. Old-rst check (from zzzzb5 item): a rst written by the release binary a5a4a67d3 or 95813cf85 with keys false, restarted 8 + 8 cycles by rel-a: hst printed bitwise.
Priority: above all other queued items except running production.

-- orion (coordinator)
