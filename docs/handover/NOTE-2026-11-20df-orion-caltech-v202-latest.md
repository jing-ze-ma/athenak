# NOTE 2026-11-20df orion -> caltech (USER): new LATEST is v20.2 (a0c19e459)

- **origin/clean-stagger2-latest = v20.2, commit a0c19e459** (docs) on f4cefb181 (code: v20.1 + the theta-wall fix + the rhse vcs_max MPI_MAX fix). It passed the gate (V202*: bitwise to V201D, rt 5/5, poison OK). Same source as B_v202.
- **Production keys gain `hydro/floor_pay = rad_column`** for all stars. On a **fresh start it must be in the deck's `<hydro>` block**; naming it only on the command line is FATAL ("not found"). Restarts accept it on the command line.
- From now on, v20.x baselines and controls use the v20.2 binary (B_v202) with floor_pay on.
- Running arms are unaffected.

-- orion (for the USER)
