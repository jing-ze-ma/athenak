# TASK 2026-10-21 deltaai -> raven: gas temperature from dumps (USER: "use gas temperature" on the RSG E DeltaAI page)

The RSG E page (claude.ai/artifact/JegJZ1mkx7nqHqz3pk2aU9; your rsg-page c542d6b0 + native-cell slices/movies) currently shows
T_rad. USER wants T_gas. hydro_w dumps hold dens + eint only; the code's <output> "temperature" derived variable is p/rho (ideal,
wrong for eos = general/table), so T_gas has to be inverted offline.
1. Do you have a Python T_gas(rho, eint) [eint per volume, as in hydro_w] that reproduces the code's general_eos = table build for
   the RSG E deck (eos_ionization true: H + He Saha; eos_h2 false; eos_metal_ionization false; eos_radiation false; xh 0.653848,
   yhe 0.326395, a_metal 16; logd -21..-1, logT 3.3..6.8, dlog 0.02; tfloor 2000 K)? If yes push it (e.g. docs/handover/rsg-page/
   tgas.py) with a check value or two. We are writing one from src/eos meanwhile and will cross-check against yours + the IC column.
2. If the code can dump its EOS table (or T) without a deck physics change, say how (an extra <output> block = output-only; I would
   ask before adding it).
3. Link 1 (3425380) ends ~20:56 CDT; the bg-check NOTE you asked for (c21e4d7a) follows right after.
