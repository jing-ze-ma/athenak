# NOTE Raven -> Caltech: He giant FIN4 binary = fin4-1014 @ 5c504e49 (if not built yet)

Addendum to TASK-2026-10-14-caltech-hegiant-fin4.md (worker HG). MG's vet_gd cost fix: fin4-1014 @ 5c504e49 is BITWISE to d505d21f
(fresh runs + old-rst continuation, A100 + MI300A) and ~3.5 % faster. Raven/Viper He giant production now run 5c504e49
(A100 b995f6e3, MI300A ec9722e8).
- If you have NOT built d505d21f yet: build **5c504e49** instead (same flags; name athena_gpu_he_star_m1_5c504e49_nofma).
- If you already built d505d21f: keep it (bitwise).
- Raven + Viper validated the FIN4 He giant (smoke + restart BITWISE on A100 and MI300A; Raven background run passes: budget
  L_in + P_esrc = 6.755e6 Lsun = L_MESA(r_in); bands identical to the FIN3 reference E0). Your section-3C background arm may be
  SKIPPED; the smoke + restart gate (3A/3B) is still required on your own build before link 1.
