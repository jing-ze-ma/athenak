#!/usr/bin/env python3
"""mk_table.py <rs_eval json> -- accuracy table of the rsx family (end state), one row per
run, plus the gate verdict.  Gates (vs the LE08 DIFFUSION solution; the IDORT transport
reference differs from it by L1pre 9-15 %, Lpre(0.01) +15-21 %, Lpre(0.1) +2.5-7 %,
T_r,s up to +2.2 % for M0 = 1.2-3, transport_vs_le08.py, so those are NOT gated tightly):
  vs LE08: L1 rho <= 2 % (t7 criterion), L1 Tg, Tr <= 3 % (M3 transport-diffusion = 2.3/2.5 %)
  |Lpre(0.1) error| <= 10 %
  |T_r at the shock error| <= 3 %
  |downstream rho, Tg, Tr| <= 1 %
  shock speed at the end |v_s| / v0 <= 1e-3 (the early relocation is the IC transient)
  vs TRANSPORT (M0 = 1.2/2/3): L1 <= 0.5 %, |Lpre(0.1)| <= 3 %, |T_r,s| <= 1 %
  mass flux deviation <= 1e-2 more than 10 cells from the shock (ringing reported apart)
  run: rc 0, no FATAL, NON-CONVERGED 0 (checked by run.sbatch output, not here)
"""
import json
import sys

G = {"L1rho": 0.02, "L1": 0.03, "Lpre": 0.10, "Trs": 0.03, "down": 0.01, "vshock": 1e-3,
     "mflux": 1e-2, "TL1": 0.005, "TLpre": 0.03, "TTrs": 0.01}


def verdict(o, lx):
    """LE08 (diffusion) gates allow the transport-vs-diffusion model difference; where the
    IDORT transport solution exists (M0 = 1.2/2/3) the tight gates are against it."""
    bad = []
    if o["L1_rho"] > G["L1rho"]:
        bad.append("L1_rho")
    for k in ("L1_Tg", "L1_Tr"):
        if o[k] > G["L1"]:
            bad.append(k)
    if abs(o["eLpre_0.1"]) > G["Lpre"]:
        bad.append("Lpre0.1")
    if abs(o["eTrs"]) > G["Trs"]:
        bad.append("Trs")
    for k in ("edown_rho", "edown_Tg", "edown_Tr"):
        if abs(o[k]) > G["down"]:
            bad.append(k)
    if abs(o.get("vshock", 0.0)) > G["vshock"]:
        bad.append("vshock")
    if o["mflux"] > G["mflux"]:
        bad.append("mflux")
    if "T_L1_Tg" in o:
        for k in ("T_L1_rho", "T_L1_Tg", "T_L1_Tr"):
            if o[k] > G["TL1"]:
                bad.append(k)
        if abs(o["T_eLpre_0.1"]) > G["TLpre"]:
            bad.append("T_Lpre0.1")
        if abs(o["T_eTrs"]) > G["TTrs"]:
            bad.append("T_Trs")
    return "PASS" if not bad else "FAIL(" + ",".join(bad) + ")"


def main():
    rows = json.load(open(sys.argv[1]))
    print("| run | kind | drift [cells] / v_s/v0 | L1 rho / Tg / Tr | L1pre | dLpre0.1 | dLpre0.01 "
          "| dTrs | Tmax code / exact / cell-avg exact | down rho / Tg | vs TRANSPORT L1 Tg, "
          "dLpre0.1, dTrs | mflux far / near | verdict |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for o in rows:
        if o.get("tag") != "end" or "error" in o:
            if "error" in o:
                print("| %s | ERROR %s |" % (o["dir"], o["error"]))
            continue
        lx = o["dx"] * o["nx1"]
        name = "/".join(o["dir"].rstrip("/").split("/")[-4:-1]) if o["dir"].endswith("out") \
            else "/".join(o["dir"].rstrip("/").split("/")[-3:])
        tr = ("-" if "T_L1_Tg" not in o else "%.4f, %+.1f%%, %+.2f%%"
              % (o["T_L1_Tg"], 100 * o["T_eLpre_0.1"], 100 * o["T_eTrs"]))
        print("| %s | %s | %+.2f / %+.1e | %.4f / %.4f / %.4f | %.3f | %+.1f%% | %+.1f%% | %+.2f%% "
              "| %.3f / %.3f / %.3f | %+.2f%% / %+.2f%% | %s | %.0e / %.0e | %s |"
              % (name, o["kind"], o["drift"], o.get("vshock", float("nan")), o["L1_rho"], o["L1_Tg"],
                 o["L1_Tr"], o["L1pre"], 100 * o["eLpre_0.1"], 100 * o["eLpre_0.01"],
                 100 * o["eTrs"], o["Tmax"], o["Tmax_ref"], o["Tmax_refcell"],
                 100 * o["edown_rho"], 100 * o["edown_Tg"], tr, o["mflux"],
                 o.get("mflux_near", float("nan")), verdict(o, lx)))


if __name__ == "__main__":
    main()
