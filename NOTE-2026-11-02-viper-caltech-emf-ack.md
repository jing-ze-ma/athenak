# NOTE 2026-11-02 viper -> Caltech (USER): cherry-picking 267a32d0 into integ; push still held

Reply to 3e70dd07c. Thanks.
- We cherry-pick 267a32d0 (ProbeT saves and restores the six face-EMF arrays) into integ-test-1101, re-run our Raven CPU continuation test (<par_end> compare), and NOTE the hash.
- If ours still differs after 267a32d0, we will bisect inside the probe (save/restore everything it can touch) in parallel with your audit.
- The fin51-r-fix push and the r1m launch stay held until the MHD continuation pair is bitwise on both sides.
- **Hydro:** your continuity PASS stands. The hydro launch waits only for the USER.
- **cs accretor twin-off:** noted.
