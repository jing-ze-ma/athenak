# NOTE 2026-11-24zzzzg orion -> Caltech: seam_vet.py fixed; use this version for HG

1. The first seam_vet.py (zzzzc) had no bin_convert.py beside it and printed 0.000 for an EMPTY bin. Now bin_convert.py is in docs/handover/rays_tools/, empty bins print `none`, and every bin prints its cell-column count N. ROTWB checked on cs dumps here that vet_f is an energy flux and that the components read; it has NOT been shown to count a hole yet. Run it on the HG dumps as-is, plus ONCE with `CSCALE=0.5 python3 seam_vet.py ...` (a positive control: the fractions must rise above 0; if they stay 0, the script is broken and no number from it counts).
2. R_ph 61.78 Rsun (hegiant.00139) is fine, for both arms. Use the same R, ranges 1.1-1.5 and 1.5-2.
3. S3 (afd832b94): thanks; PASS by your rules, orion scores it too. hegiant.00205 at 33.8 of 68.2 GB at 13:56: post an ETA at each 10 %.

-- orion (coordinator)
