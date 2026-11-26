# NOTE 2026-11-27zx orion -> Caltech: scpr4, rerun 5b only (curved J output was the wrong variable)

Thanks for 27zw; everything else is done. 5b gave J = 0 because `vet_vsc` holds the Cartesian rays
only. On the curved path, the rays' J is `vet_tj` in `variable = vet_thsw` (-1 below the ray cut).
`sph_err.py` now masks those cells (pull scpr-1127).

Please rerun 5b only, 1 GPU, same binary (rel6 suite 8cc98fb3):
`-i our_sph_sp.athinput output1/variable=vet_thsw`, then
`SIG=50 sph_err.py <last bin> vet_tj sp`. Post the printed lines. If every cell is -1, post
`grep -i "vet_gd\|cut\|thsw" log | head -20` instead and stop. Delete the bin. Reply caltech -> orion "scpr4".

-- orion (SCPR, for the coordinator)
