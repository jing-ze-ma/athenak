---
name: fastchem-vs-general-eos
description: "MEASURED 2026-08-22: the FastChem CE table buys the general EOS nothing -- mu agrees to <=2%, x_e to <=25% wherever it matters, and the table cannot supply internal energy at all. DECIDED against; do not re-propose."
metadata:
  node_type: memory
  type: project
---

Asked whether `data/exo_fms_ck/CE_tables/FastChem_ck_1x_int.txt`, already fetched for the
correlated-k continuum, could improve the general EOS. **Answer: no.** Measured, not
argued -- `eos_composition.hpp` is host-only plain C++ (no Kokkos), so it compiles into a
standalone driver against `-I src` and can be evaluated on the FastChem grid directly.
Invert `p = rho * p_spec(rho,T)` for rho by bisection at each (T,p) node.

**What the table actually holds:** 241 T (100-6100 K, 25 K spacing) x 34 p (1e-8..1000
bar), six columns per record -- mu, then VMR of H2, He, H, e-, H-. FastChem tracks 33
species internally but **only these five are written out, and there is no internal
energy**. So it cannot supply `e(rho,T)`, which is the whole job of the general EOS: the
H2 dissociation reservoir is 1.58e12 erg/g and mu alone does not reconstruct it.

**Comparison over the dhj domain (4e-6..250 bar), AthenaK / FastChem:**

| quantity | agreement |
|---|---|
| mu, everywhere | median 1.005, worst 0.979 (6100 K, 214 bar) to 1.006 |
| x_e, 1500-2500 K | median 1.05 |
| x_e, 2500-3500 K | median 1.07 |
| x_e, 3500-4800 K | median 0.93 |
| x_e, > 4800 K | median 0.94, worst 0.83 |
| x_e, 1000-1200 K | median 1.45 |
| x_e, 800-1000 K | median 290 |
| x_e, < 800 K | 1e5 - 1e10 |

**The low-T blowup is irrelevant, and that is the decisive point.** `max_eta` caps the
resistivity and therefore FLOORS x_e at `230 sqrt(T)/max_eta`. Below 1500 K, 89 % of
domain cells have BOTH codes under that floor at max_eta = 1e14 (95 % at 1e13), and the
fraction of cells where the two disagree by more than 2x AND either is above the floor is
**0.0 %**. Swapping in FastChem's x_e would not move the resistivity anywhere.

**Coverage kills it independently:** the table stops at 6100 K, but the dhj deep interior
is 5100-12000 K with a T 90th percentile of 10113 K, and 40.6 % of cells exceed 4000 K.
The hot gas where the resistivity is largest is off the top of the table. The EOS's own
table spans log10 T 1.5-6. It is also a (T,p) table being asked to serve a (rho,T) EOS.

**The real value of this exercise is the VALIDATION, and it is worth citing:** the
general EOS's mean molecular weight agrees with an independent 33-species equilibrium
chemistry code to <= 2 % over the whole domain, and its electron fraction -- Saha with
singly-ionizing metal donors and per-species rainout -- to <= 25 % over 1500-6100 K.
That is a strong independent check on [[xe-resistivity-long-runs]]'s premise.

**Written up in `docs/general_eos.md` under "Checked against equilibrium chemistry"
(9a8fe0bd)**, cross-linked from `docs/correlated_k_rt.md`. The same commit fixed the
validity section, which still listed metal condensation as not included even though the
section above it documents how to switch it on.

See [[correlated-k-design]] for where the table came from.
