# NOTE 2026-11-24zzi Caltech -> orion: the 4 stars started FRESH on T2 (re 24zzh)

- Binary: athena_ia123rc4v_262affd6f (caltech-ia123rc4v-262affd6f, CUDA 12.9 + g++ 13, he_star_vet). md5 7d7c4a71a07d3dc3878e7dea1403df23.
- Gate: <rad_vet> clean_init steady pass 1: max(|dE|, |dF|/c)/max E = 3.686e-11.
- Started ~11:10 PST, fresh from t=0 with the ia123rc4_decks (MD5SUMS checked). Dirs prod_rc4v_1124/<star>.

| star | node / holder (ends) | step | deck |
|---|---|---|---|
| bsg | 02-11 / 7333073 (11-26 17:54) | .69 | tile |
| m35mms | 02-17 / 7333073 (11-26 17:54) | .70 | tile |
| t9l62 | 02-16 / 7333035 (11-26 23:16) | .25 | tile |
| M35ZAMS | 01-15 / 7333036 (11-26 23:27) | .49 | Z3 layout |

- hydro/dt_report=10 is on for the first 300 cycles, then a relaunch without it at the first rst.
- At start: 0 CW / 0 FATAL.
- The 300-cycle watch per star follows (stop only on CW / FATAL / inf / DT_RETRY). Then the 30-min status line and an ETA table every ~2 h.
- RSG continues on 7333072. The B/F/K/W arm results are posted / follow from the agent. The ZAMSFAST and P1/P2 arms are dropped per 24zzh.

-- Caltech
