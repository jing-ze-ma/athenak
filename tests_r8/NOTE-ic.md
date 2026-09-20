# The ic and its generator live OUTSIDE the repository

`inputs/hydro/he4_presn_cs.athinput` sets

    problem/ic_profile = /viper/u2/jinma/ATHENAK/bench/hestar_presn/ic_he4_presn_sph.txt

and that file, its generator `make_ic_sph.py` and the stellar column
`column_he4_presn_sph.txt` are in `bench/hestar_presn/`, which is not a git repository.
The generator as it stands after this round is committed here as `make_ic_sph.py` (a COPY,
for the record -- edit the one in `bench/hestar_presn/`, which is what produced the file
the input points at).  The ic itself is 800 kB of numbers and is not committed; regenerate
it with

    python3 /viper/u2/jinma/ATHENAK/bench/hestar_presn/make_ic_sph.py

The pre-gate ic is preserved nowhere: it is `T_HI = 0.0` in that script.
