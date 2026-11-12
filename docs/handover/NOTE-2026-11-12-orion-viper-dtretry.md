# NOTE 2026-11-12 Orion -> viper: dtretry CPU gate (HG small cs deck)

Builds (PROBLEM=he_star_m1, Release, MPI, -ffp-contract=off, kokkos from the Orion src_vfix copy):
- base 84cd20888: athena md5 004d17ab403edbd3c34ec62cbc4750c6
- fix 3ec401877: athena md5 e65846c6c5bcd232bcbefc4836ea7415

Tables are from docs/handover/hgpos-1109/files. The MD5SUMS file is not on the branch, so the tables were not checksum-checked; they are the files
already used for the hgpos/hgfix gates. Ran the posted run.sh unchanged: orion003, 7 arms x 3 ranks, wall 348 s.

## Result
| check | expected | CPU result |
|---|---|---|
| BN vs CN (keys absent) | bitwise | **IDENTICAL** (00002-00004) |
| CN vs CP (forced retry, views=all, audit, poison=1) | bitwise | **IDENTICAL**; 0 FATAL, 0 POS_TRIGGER |
| CTN vs CTP (TIME keys, same) | bitwise | **IDENTICAL**; 0 FATAL |
| CN vs CF (forced retry, views minimal = production) | bitwise | **DIFFER** from 00003 on (00002 IDENTICAL): 5,726,605 of ~175M words, maxrel 1 |
| CTN vs CTF (TIME keys, views minimal) | bitwise | **DIFFER** from 00003 on: 5,726,539 words, maxrel 1 |

**So on CPU the fix is clean with views=all, but not with the production minimal list.** The retry at cycle 3 changes the run in both key sets.
This differs from viper's GPU result (forced retry == straight with minimal views).
- **The 51 unchecked Views are not the cause.** poison=1 NaNs the "51 changed Views outside the minimal list", and CP/CTP stay bitwise with views=all. With poison on top of views=all, those 51 are rewritten before they are read.
- **Something only the full checkpoint restores is read across the retry on CPU.** The audit lists 72 (prod) / 75 (TIME: + fax1..3) checkpointed Views changed by the failed attempt; the difference set (all minus minimal) is the candidate list.
- Note: CPU vs GPU may also differ in which scratch Views a kernel writes. The CPU (OpenMP/serial) branches fill some Views that the GPU path computes on the fly, e.g. per-cell caches under od_cache / krylov_fuse being off.

## gate.out (whole)
```
BN rc=0 fatal=0 postrig=0  
CN rc=0 fatal=0 postrig=0  
CTN rc=0 fatal=0 postrig=0  
CTF rc=0 fatal=0 postrig=0  
CF rc=0 fatal=0 postrig=0  
CTP rc=0 fatal=0 postrig=0 <time> dt_retry poison: 51 changed Views outside the minimal 
CP rc=0 fatal=0 postrig=0 <time> dt_retry poison: 51 changed Views outside the minimal 
== BN vs CN
hegiant.00002.rst IDENTICAL
hegiant.00003.rst IDENTICAL
hegiant.00004.rst IDENTICAL
== CN vs CF
hegiant.00002.rst IDENTICAL
hegiant.00003.rst DIFFER len 1403580406/1403580406 ndiff 5726605 maxrel 1.000e+00 (8-byte words, header-aligned approx)
hegiant.00004.rst DIFFER len 1403580406/1403580406 ndiff 5726605 maxrel 1.000e+00 (8-byte words, header-aligned approx)
== CN vs CP
hegiant.00002.rst IDENTICAL
hegiant.00003.rst IDENTICAL
hegiant.00004.rst IDENTICAL
== CTN vs CTF
hegiant.00002.rst IDENTICAL
hegiant.00003.rst DIFFER len 1502548758/1502548758 ndiff 5726539 maxrel 1.000e+00 (8-byte words, header-aligned approx)
hegiant.00004.rst DIFFER len 1502548758/1502548758 ndiff 5726539 maxrel 1.000e+00 (8-byte words, header-aligned approx)
== CTN vs CTP
hegiant.00002.rst IDENTICAL
hegiant.00003.rst IDENTICAL
hegiant.00004.rst IDENTICAL
wall 348 s
```

## grep -E "DT_RETRY|dt_retry (poison|audit)" runs/<arm>/run.log | cut -c1-400 (the DIFFER arms CF, CTF)
```
--- CF
DT_RETRY cycle=3 t=4.387500e+02 dt 1.462500e+02 -> 1.462500e+02 reason=forced resid=0.000000e+00 cut 1
--- CTF
DT_RETRY cycle=3 t=4.387500e+02 dt 1.462500e+02 -> 1.462500e+02 reason=forced resid=0.000000e+00 cut 1
```

## Full audit lines (CP, CTP) for the candidate list
```
--- CP
<time> dt_retry poison: 51 changed Views outside the minimal list set to NaN
<time> dt_retry audit cycle 3: 72 checkpointed Views changed by the failed attempt: hyd.u0(2.749120e+01) hyd.w0(2.749120e+01) hyd.wder(1.099648e+01) hyd.wtemp(5.498240e+00) hyd.u1(2.749120e+01) hyd.uflx.x1f(2.749120e+01) hyd.uflx.x2f(2.749120e+01) hyd.uflx.x3f(2.749120e+01) hyd.utest(2.749120e+01) hyd.wbq0(8.247360e+01) hyd.sbf_s(2.749120e+01) m1.esrc(5.498240e+00) m1.pos_cnt_d(8.000000e-05) m1.cf_r(5.498240e+00) m1.f0x1(5.505984e+00) m1.f0x1n(5.505984e+00) m1.iw(5.058381e+02) m1.ktn(5.498240e+00) m1.ifw(2.202394e+01) m1.vgd_hr(5.498240e+01) m1.ifw2(2.299264e+01) m1.ifw3(2.299264e+01) m1.gcz1(5.505984e+00) m1.gcz2(5.748160e+00) m1.gcz3(5.748160e+00) m1.csk1(1.649472e+01) m1.csk2(1.649472e+01) m1.csinc(1.649472e+01) m1.f0x2(5.748160e+00) m1.f0x3(5.748160e+00) m1.f0x2n(5.748160e+00) m1.f0x3n(5.748160e+00) m1.thw(7.697536e+01) m1.thq(2.199296e+01) m1.krw(5.498240e+00) m1.ost(1.044666e+02) m1.lay_part(2.097216e+00) m1.lay_part_h(6.400000e-05) m1.ipred(2.749120e+01) m1.t2gsv_h(2.749120e+01) m1.dbrow(8.797184e+01) m1.t2k1(4.398592e+01) m1.t2k2(4.398592e+01) m1.t2inc(4.398592e+01) m1.t2f1(5.505984e+00) m1.t2f2(5.748160e+00) m1.t2f3(5.748160e+00) m1.ipred2(2.749120e+01) m1.vet_now(1.099648e+01) m1.vet_opac(2.199296e+01) m1.x2_tt(5.498240e+01) m1.x2_hr(5.498240e+01) m1.x2_q(7.744000e-03) m1.thx2(5.748160e+00) m1.thx3(5.748160e+00) m1.vmw(6.597888e+01) m1.d3w(1.649472e+01) m1.amw(3.298944e+01) m1.u0(2.199296e+01) m1.u1(2.199296e+01) m1.opac(2.199296e+01) m1.tau_ten(5.498240e+01) m1.vcol_q(7.744000e-03) m1.vlat_cs(1.099648e+01) m1.vlat_lx(3.848768e+01) m1.vra_a0(5.498240e+01) m1.vra_a1(5.498240e+01) m1.vra_c(7.744000e-03) m1.vgd_rg_store(6.935560e+00) m1.vgd_cs(2.944512e+01) m1.vgd_cssb(3.225600e-01) m1.vgd_csrb(3.225600e-01)
--- CTP
<time> dt_retry poison: 51 changed Views outside the minimal list set to NaN
<time> dt_retry audit cycle 3: 75 checkpointed Views changed by the failed attempt: hyd.u0(2.749120e+01) hyd.w0(2.749120e+01) hyd.wder(1.099648e+01) hyd.wtemp(5.498240e+00) hyd.u1(2.749120e+01) hyd.uflx.x1f(2.749120e+01) hyd.uflx.x2f(2.749120e+01) hyd.uflx.x3f(2.749120e+01) hyd.utest(2.749120e+01) hyd.wbq0(8.247360e+01) hyd.sbf_s(2.749120e+01) m1.esrc(5.498240e+00) m1.pos_cnt_d(8.000000e-05) m1.cf_r(5.498240e+00) m1.f0x1(5.505984e+00) m1.f0x1n(5.505984e+00) m1.iw(5.058381e+02) m1.ktn(5.498240e+00) m1.fax1(5.505984e+00) m1.fax2(5.748160e+00) m1.fax3(5.748160e+00) m1.ifw(2.202394e+01) m1.vgd_hr(5.498240e+01) m1.ifw2(2.299264e+01) m1.ifw3(2.299264e+01) m1.gcz1(5.505984e+00) m1.gcz2(5.748160e+00) m1.gcz3(5.748160e+00) m1.csk1(1.649472e+01) m1.csk2(1.649472e+01) m1.csinc(1.649472e+01) m1.f0x2(5.748160e+00) m1.f0x3(5.748160e+00) m1.f0x2n(5.748160e+00) m1.f0x3n(5.748160e+00) m1.thw(7.697536e+01) m1.thq(2.199296e+01) m1.krw(5.498240e+00) m1.ost(1.044666e+02) m1.lay_part(2.097216e+00) m1.lay_part_h(6.400000e-05) m1.ipred(2.749120e+01) m1.t2gsv_h(2.749120e+01) m1.dbrow(8.797184e+01) m1.t2k1(4.398592e+01) m1.t2k2(4.398592e+01) m1.t2inc(4.398592e+01) m1.t2f1(5.505984e+00) m1.t2f2(5.748160e+00) m1.t2f3(5.748160e+00) m1.ipred2(2.749120e+01) m1.vet_now(1.099648e+01) m1.vet_opac(2.199296e+01) m1.x2_tt(5.498240e+01) m1.x2_hr(5.498240e+01) m1.x2_q(7.744000e-03) m1.thx2(5.748160e+00) m1.thx3(5.748160e+00) m1.vmw(6.597888e+01) m1.d3w(1.649472e+01) m1.amw(3.298944e+01) m1.u0(2.199296e+01) m1.u1(2.199296e+01) m1.opac(2.199296e+01) m1.tau_ten(5.498240e+01) m1.vcol_q(7.744000e-03) m1.vlat_cs(1.099648e+01) m1.vlat_lx(3.848768e+01) m1.vra_a0(5.498240e+01) m1.vra_a1(5.498240e+01) m1.vra_c(7.744000e-03) m1.vgd_rg_store(6.935560e+00) m1.vgd_cs(2.944512e+01) m1.vgd_cssb(3.225600e-01) m1.vgd_csrb(3.225600e-01)
```

Run trees (restart files, ~1.4 GB each): Orion /orion/ptmp/jinma/dtretry_1112/dtretry_1111/runs; kept for a follow-up bisection (e.g. views = minimal + one Group) unless you say delete.
