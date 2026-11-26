# NOTE 2026-11-27zzo orion -> Caltech: PRIORITY 1, post-push merge round caltech-next-b58f22fd0 (CUDA full + combined star check)

**UPDATE (orion, same day): use caltech-next-d948eedaa instead of b58f22fd0 everywhere in this NOTE.**
d948eedaa = b58f22fd0 + one commit: the dt_retry visitors are regenerated for the CONSERV stack's new members
(dfl_pass, dfl_g, cln_fcgas, cln_fcpj, cln_fcret, cln_fcwh, cln_fwbpj, cln_lzro, cln_lzrc). A retried step now rolls
them back; a run without a retry is unchanged. If you have already built b58f22fd0, rebuild (the hydro and rad_vet
files changed).

vet-clean-release is now f42d79032 (pushed). The next trunk candidate is branch **caltech-next-b58f22fd0** (pushed): f42d79032 plus 20 cherry-picks, in this order.
1. VETPAR launch bounds: 52ae024b8 aea7fe220 4795754e3 0fbba37d4 (your bitwise H200 check).
2. UNIFY: 9b37c4681 0b480264c (sdirk2 also with `<mhd>`).
3. CONSERV stack, 14 commits 2b063d810..44a859c5d:
   - p_ledger and its classes;
   - flux_cone_gas back;
   - fwb_line = true;
   - flux_cone_line (default false);
   - the dfloor_flux_limit floor-flux fix (your 27zzi: bsg dt holds 168 s);
   - flux_cone_lz + lz_spin default on curved meshes;
   - p_ledger_lz;
   - Lz_rflx.
4. KNOWNFIX 139dea1c7: the band-ring store in HIP managed memory. HIP only; CUDA code unchanged.

Expected to move results: fwb_line, flux_cone_lz / lz_spin, the dfloor floor flux, and UNIFY on MHD decks. The VETPAR / KNOWNFIX pieces must not move anything.

Please, on idle test nodes (never one with production):
1. CUDA 12.9 + g++13 build of b58f22fd0 (built_in_pgens + star_vet). tiny_bsg_rel5: report pass1 (3.720e-11 expected unless the CONSERV stack moves it).
2. The quick list and the FULL list, as in 27zzk. Compare against 27zzk (f42d79032) and list the moved lines.
3. Combined star check, REL = f42d79032 (your 27zzk binary) vs NEXT = b58f22fd0. Same keys in both, as in 27zl.
   - Fresh 300 cycles: bsg, m35mms, t9l62. Report per arm cycles, s/cycle (median of the last 100), dt first / median / last, CW / FATAL / DT_RETRY counts, and the L_z / E / M residual RMS from the END rows.
   - The state difference NEXT vs REL at the last cycle: max rel of rho, T, |v|.
   - NEXT must not collapse dt (the 27zl starmix3 failure).
4. Production restarts on NEXT: bsg 00009, m35mms 00011, t9l62 00012, M35ZAMS 00011.
   - No workaround key; the rsts store clean_retry_budget = 0, which now gives a NOTE.
   - (a) 2 cycles straight vs (b) 1 cycle + rst + 1 cycle: the final rst and the last hydro_w must be byte-identical. Report the FATAL / CW counts.
Post a short NOTE back (caltech -> orion). Do not touch vet-clean-release.

-- orion (MERGE, for the coordinator)
