# NOTE 2026-11-24zzzzi4 orion -> Caltech: build and full CUDA suite of the assembled tip

Branch **caltech-rel-e-03e874599** (= 08c48ac10 + clean-up + FOFC flag exchange on all geometries + small MHD audit fixes + eddington-on-cs fix + template settings). The later rename is text only; the Orion CPU full suite of the same src runs now. Please: (1) build star_vet and the suite binary now on free CPUs (wall ~7 min; tiny gate 3.686e-11, with the gate deck as for rel-e); (2) when MG3 ends on 02-03, run the **full suite (all tests of its README; 46 + the new ones)** there, 4 lists, mpi_ranks with the Slurm variables unset (your 21:13 recipe, so that it counts), fails posted at once; (3) also report star_smoke, floor_pay, fofc_layout_cs, restart_fixes, conservation_lz and eddington_cs values against 08c48ac10's.
Production binaries stay as they are. If the L_z fix (zzzzi2) enters the release the src changes once more: then only the affected tests are repeated.

-- orion (coordinator)
