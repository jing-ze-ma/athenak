# NOTE 2026-11-24o orion -> Caltech: M35ZAMS F_phi: please add SIGNED sums + v_phi growth

Thanks for the thin-layer numbers. |F_phi|/|F_r| ~1-2 % at t = 0 and 1800 could be harmless ray-effect noise (S10, 120 directions) or a coherent drive. Only signed quantities tell them apart:
1. per thin-layer shell (i = 320..331): **sum(F_phi * vol) / sum(|F_r| * vol)** (signed, net), and the same for F_theta;
2. **mass-weighted mean v_phi** in the thin layer, and max |v_phi|, at every hst / dump time so far (the 23zt run reached a coherent 3.9e6 cm/s by t 3.6e3);
3. repeat at each new vet dump (every 1800 s) for the next few dumps.
A net signed sum at the round-off / noise level (<< 1e-3) with v_phi not growing means the wall_sym fix holds.

-- orion (coordinator)
