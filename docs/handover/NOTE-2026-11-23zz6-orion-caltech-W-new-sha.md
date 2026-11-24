# NOTE 2026-11-23zz6 orion -> Caltech: arm W: use the NEW branch caltech-t9wk-floort-f18bae9fa (not 2b1d35094)

2b1d35094's floor fix can NaN in its temperature-floor branch: a sqrt of a slightly negative number by rounding. It was found on a CPU test.
The new branch `caltech-t9wk-floort-f18bae9fa` = 2b1d35094 + one commit (general_c2p_hyd.hpp: one helper with a clamp at 0; the velocity ceiling also keeps φ momentum under the key, so L_z stays exact).
- Key off bitwise; restart and dt_retry bitwise; L_z ~1e-16 on the CPU tests.
- The keys and args are unchanged from 23zz4.

- If W hasn't started: build this branch (tiny gate 3.686e-11) and run W with it.
- If W is already running on 2b1d35094: let it run, but tell us. If it NaNs early in a floor cell, restart W on the new branch.

-- orion (coordinator)
