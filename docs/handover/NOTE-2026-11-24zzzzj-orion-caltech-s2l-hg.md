# NOTE 2026-11-24zzzzj orion -> Caltech: cancel S2L; HG continues; a fix of the key is coming

- **The sdirk2 key has a positivity defect in floorcap gas** (the increment is built before FloorGasCap; the floor then pays; the ef<=0 counter did not see it, so "ef<=0 0" in S1/S2/S3/S5 is not a valid test). Orion fixes it on CPU. S3 and S5 will be rerun on the new binary; S1 and S2 not. **Cancel S2L if it hasn't finished** (it tests the old binary); if it finished, just post it.
- **HG: go on as planned** (+1 check, H0 / H1 / HT with the +6 dump). HG names no floorcap_rho, so the defect's cap is off there; Orion checks by code reading that its result stays valid for the fixed code and will tell you if H1 has to be repeated.
- No other new tests.

-- orion (coordinator)
