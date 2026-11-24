# NOTE 2026-11-23zz7 orion -> Caltech: speed up the test queue (USER: progress too slow)

To get answers sooner:
1. **Stop WA now** (WB contains A plus the face clamp, so WA adds nothing). **Start W on WA's node immediately**, with branch caltech-t9wk-floort-f18bae9fa (23zz6). W is the decisive t9l62 test.
2. Keep WB and SB running to cycle 8050.
3. On the next free node: the m35mms nsys profile (23zz3, ~10 min), then the t9l62 DTREP arms T0 / T0d (23zz5, ~30 min each). Run two at a time on one node if memory allows (m35mms and the DTREP arms use < 40 GB / GPU).
4. Post a one-line status per arm every 30 min (cycle, t, dt, CW, fallbacks, DT_RETRY), so we can react without waiting for the end.

m35mms production on IA11: unchanged.

-- orion (coordinator)
