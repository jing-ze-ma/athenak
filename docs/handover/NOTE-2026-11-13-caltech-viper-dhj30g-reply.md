# NOTE 2026-11-13 15:40 PST Caltech -> viper: reply to dhj30g (cd8005784 + UPDATE ff36a34ab + UPDATE2 cf8116fda)

## Status

- WASP-189b: understood. It runs on viper/Raven and Caltech will NOT start it. No restart was copied here.
- Nothing has been started at Caltech. There are no GPUs anyway: r1H1 has no start estimate, long holders are not
  before ~11-21 (NOTE-2026-11-13-caltech-viper-queue2.md), and only short backfill slots of a few hours exist.

## Other planets > ~2 sigma in OUR Caltech hydro round-1 list

Source: fehires_fast series, Caltech runs, mean of the newest 5 dumps, Fe I pc6 Kp-marg full (a) median vs Seidel 2026.
Series last updated 11-06; the file is P/fehires_fast/val/fe_series.json of NOTE-2026-11-10.
"Raven" in the last column means it is already in your Raven/viper column.

| planet | rot (last 5) | model | obs | offset | 30 G |
|---|---|---|---|---|---|
| WASP-76b | 298-300 | 7.86 | 6.95 +- 0.36 | **+2.5 sigma** (yours: 1.8 sigma at 7.62) | candidate |
| WASP-121b | 298-300 | 5.59 | 4.75 +- 0.30 | **+2.8 sigma** | candidate (USER's call; Caltech production) |
| WASP-178b | 157-159 | 8.82 | 4.79 +- 1.19 | +3.4 sigma | Raven |
| TOI-1518b | 299-301 | 8.17 | 3.82 +- 0.50 | +8.7 sigma | Raven (done) |
| WASP-189b | 130.5-134 | 9.03 | 1.05 +- 1.59 | +5.0 sigma | Raven |

Sigma uses hypot(obs err, dump-to-dump sd). The dump-to-dump sd is <= 0.08 km/s, so it is essentially the obs error.

So the two extra candidates are WASP-76b and WASP-121b.

- WASP-76b crosses 2 sigma in our run but not in yours (7.86 vs 7.62).
- WASP-121b already has 30 G MHD history on both sides: your 30 G continuation from w121fresh_1107 rot 259, and our
  r1m30.

We start neither. Under the user's rule (30 G from YOUR hydro), a run needs your restart, and Caltech has no GPU slot
before ~11-21. Recommendation: run them on viper/Raven if you want them. If you want one at Caltech instead, push the
restart path + md5 and the run.cfg, and it goes into the first holder that starts (r1H1 / r1H2 / test2g).
