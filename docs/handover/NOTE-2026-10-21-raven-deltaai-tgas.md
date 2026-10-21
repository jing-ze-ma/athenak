# NOTE 2026-10-21 raven -> deltaai: T_gas for the RSG E page (answers TASK 6813d9b1)

1. Yes: docs/handover/rsg-page/tgas.py inverts the CODE'S OWN general-EOS table (no re-implementation of Saha needed).
2. The code dumps its EOS table at setup with <hydro>/eos_table_dump (output-only). Do NOT add it to the production deck: run
   docs/handover/rsg-page/eos_dump_rsgE.athinput (1-cell built-in shock_tube, same EOS keys as RSG E: table, ionization on, h2/metal/
   radiation off, X 0.653848, Y 0.326395, a_metal 16, logd -21..-1, logT 3.3..6.8, dlog 0.02) with a binary built WITHOUT a user
   pgen (plain build; the he build cannot run shock_tube). Raven's dump of the same keys: md5 fef11ed03e2ae70ebf408c3c882b558a
   (24 MB, 876 x 1001 nodes; EOS code unchanged since) — compare yours.
3. Check values (round trip through tgas.py, cgs; the RSG deck units are cgs so hydro_w dens/eint go in directly):
   rho 1e-9, T 10000 K -> eint 6.443516e+03 erg/cm3; rho 1e-12, T 5012 K -> 5.177552e-01; rho 1e-6, T 316228 K -> 7.558848e+07;
   rho 1e-15, T 3162 K -> 2.886061e-04 (all back to T within 0.2 K). Apply tfloor 2000 K yourself if wanted.
