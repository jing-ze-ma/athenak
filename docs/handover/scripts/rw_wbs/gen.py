"""wbsgate-1112: radwave lists, 18 cells (P 0.0333/0.333/3.33 x tau 0.01..1000) x N16..256 (nt 2N) x edd/vetf,
arms PROD (arm=PROD19), PRODW (+wbsolve), FINAL (arm=TIME + zold/cmem_rst), FINALW (+wbsolve).
all_*.txt = the Orion lists; vip_*.txt = the viper subset (all edd; vetf without tau 0.01 N>=128 and tau 0.1 N256)."""
TOL = "amp=1e-06 rm_implicit_tol=1.0e-12 rm_implicit_lin_tol=1.0e-13 rm_implicit_tol_rel=0 rm_implicit_lin_cnorm=1.0e-13"
ARM = {"PROD": "arm=PROD19 ts=hesdirk2 vimp=true " + TOL,
       "FINAL": "arm=TIME ts=hesdirk2 vimp=true " + TOL + " rm_implicit_physflux_zold=fc rm_implicit_hr_cmem_rst=true"}
W = " rm_implicit_cs2_wbsolve=true"
PS = ["0.03333333333333333", "0.3333333333333333", "3.3333333333333335"]
TS = ["0.01", "0.1", "1", "10", "100", "1000"]
for a in ("PROD", "FINAL"):
    for w in ("", "W"):
        fa = open("all_%s%s.txt" % (a, w), "w"); fv = open("vip_%s%s.txt" % (a, w), "w")
        for p in PS:
            for t in TS:
                for n in (16, 32, 64, 128, 256):
                    for cl in ("edd", "vetf"):
                        ex = "ny=4 ref=edd" if cl == "edd" else "ref=vetqs_d"
                        ln = "%s %s x1 %d %s impl %d %s %s%s\n" % (p, t, n, cl, 2 * n, ARM[a], ex, W if w else "")
                        fa.write(ln)
                        heavy = cl == "vetf" and ((t == "0.01" and n >= 128) or (t == "0.1" and n == 256))
                        if not heavy:
                            fv.write(ln)
