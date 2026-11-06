# NOTE 2026-11-06 viper -> Caltech (USER): 1-D verdict, no collapse at his settings

**Result: 1-D at his settings does not reproduce the collapse. No inflow or KE jump near 18 d, or at any time to 30 d, on the uniform grid or his radial grid, seeded or not.** No new runs since 11-06m; this is the bottom line. Details and files are in NOTE-2026-11-06k and NOTE-2026-11-06m. All comparisons are against his 3-D Fig. 2 read by eye.

Settings: his code, units, tables, GM 260.795, cfl_rad 0.025, nzeta 5. Arms: `his1d_30d` (uniform 384, seeded), `his1d_g` (his grid, unseeded), `his1d_gs` (his grid, seeded).

| Feature | His 3-D, by eye | 1-D, all three arms |
|---|---|---|
| KE near 18 d | rises to ~5e43 erg | at most 1.8e40 at 16-20 d |
| Peak KE | ~5e43 at ~18 d | 3.3e40 to 1.5e41 at 3.9-4.0 d |
| Step in energies at 17-20 d | yes | no, steady loss |
| R(1e-10 g/cc) | shrinks after the event | stays ~51-53 Rsun to 30 d |

**What does match**
- No wall film and no Jacobi cap hits.
- Oscillation at 3-6 d with similar timing; on his grid it reaches the bottom of his 1-5e41 range. It fades after 6-8 d, while his keeps growing.
- `his1d_g` energy loss follows his to ~16 d. By 30 d it has also drifted close to his final values, but steadily rather than in a step, so this may be coincidence.

**Limits**
- A slow wall drift remains. First-face N^2 at 30 d is -68 (`his1d_g`) and -157 (`his1d_30d`), in 1e-10 s^-2.
- The older cfl_rad 1 arms (wall film, A_tc4 NaN) are not his setup and are not his collapse.
- GM of the published run is unconfirmed. At 1 d his KE (~7e40) is nearer our 20 Msun arm (5.2e40) than 1-D at his deck GM (5.2e39).

**Interpretation, not proven:** the inflow needs 3-D or something else absent from 1-D. It could also mean his deck differs from the published run. The Orion 3-D run is the test. It needs ~3.3 TB at his dump cadence to reach 18 d; its disk guard currently stops it at about 10-11 d.
