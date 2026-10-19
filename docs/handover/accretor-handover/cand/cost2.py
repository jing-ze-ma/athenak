#!/usr/bin/env python3
"""ESTIMATE: same scaling as cost.py (RY Per stage-1 basis) for the massive-gainer candidates."""
import importlib.util, sys
spec = importlib.util.spec_from_file_location("c", "/viper/ptmp2/jinma/accretor_1006/cand/cost.py")
src = open(spec.origin).read().split("sysl = [")[0]
exec(src)
ref = metric(6.8636, 6.24, 1.69, 30.3, 4.06, 6250)
rows = [  # name P Macc Mdon a Racc Tdon
 ("phiPer_Md3.5", 18.1, 7.5, 3.5, 64.6, 3.35, 8000),
 ("phiPer_Md2.5", 34.2, 8.5, 2.5, 98.5, 3.61, 8000),
 ("phiPer_Md1.8", 72.2, 9.2, 1.8, 162.2, 3.79, 8000),
 ("BY_Cru_R4.5", 106.4, 9.1, 1.7, 208.9, 4.5, 11000),
 ("W_Cru", 198.5, 8.2, 1.6, 306.0, 4.0, 5500),
 ("HD170582", 16.871, 9.0, 1.9, 61.2, 5.5, 8000),
 ("Plaskett_now", 14.396, 40.8, 5.9, 89.8, 10.4, 31000)]
print(f"{'name':13s} {'nr':>5s} {'cyc/orb':>8s} {'min/orb':>8s} {'strm_cells':>10s} {'coarsen f':>9s} {'min/orb@f':>9s}")
for n, P, Ma, Md, a, R, T in rows:
    nr, X, wc, dL1 = metric(P, Ma, Md, a, R, T)
    cyc = 1.29e5*X/ref[1]; cost = 23*(nr/448)*(cyc/1.29e5)
    f = max(1.0, wc/17.4)
    print(f"{n:13s} {nr:5.0f} {cyc:8.2e} {cost:8.0f} {wc:10.1f} {f:9.2f} {cost/f**3:9.0f}")
