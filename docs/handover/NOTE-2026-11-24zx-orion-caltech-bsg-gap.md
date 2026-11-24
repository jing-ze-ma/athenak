# NOTE 2026-11-24zx orion -> Caltech: bsg advanced only 2800 cycles between 08:10 and 09:51: was it stopped?

08:10 table: bsg c48200, t 3.720e6, 0.743 s/cyc. 09:51 status: c51000, t 4.181e6, 0.788 s/cyc. 2800 cycles x 0.79 s is ~37 min, but ~100 min passed. Please check (sacct / the chunk logs) for a stop, an idle gap, a stale STOP, or rst-write time, and post the elapsed= of each chunk since 08:00. If it's an idle gap, apply the same fix as t9l62 (a watch alert on any short step end). The others look right: m35mms 42 %, t9l62 21 % with dt recovering to 69, M35ZAMS 23 % at 0.28 s/cyc.

-- orion (coordinator)
