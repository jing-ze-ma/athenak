---
name: eostest-input-inventory
description: "What the ideal/general comparison inputs in /orion/u/jinma/ATHENAK/eostest/ cover, and how to run them"
metadata: 
  node_type: memory
  type: reference
  originSessionId: 85cd606f-82a3-4a02-a9d2-5e410a10595a
  modified: 2026-08-13T12:59:22.797Z
---

`/orion/u/jinma/ATHENAK/eostest/` holds the hand-run ideal-vs-general comparisons for the
general EOS work (83 `.athinput` files as of 2026-08-09). NOT in the repo and NOT in
[[run-directory-untouchable]] — safe to add to. Each is an `X_ideal` / `X_gen` pair; the
`_gen` file differs only by a `<units>` block, and is run with `<fluid>/eos=general`:

```bash
export OMP_NUM_THREADS=1 OMP_PROC_BIND=false     # login node oversubscribes otherwise
athena -i X_ideal.athinput -d OUT_i
athena -i X_gen.athinput   -d OUT_g hydro/eos=general
```

| prefix | pgen (`-DPROBLEM=`) | covers |
|---|---|---|
| `wb_iso_*`, `wbs_plain_*`, `wbs_perturb_*` | `hse_atm` | well-balanced, all `wb_option`s, static + perturbed |
| `cool_*`, `coolm_*` | built-in `shock_tube` | ISM cooling, hydro / MHD |
| `smr_*`, `smrm_*`, `amr_*` | built-in | static + adaptive refinement |
| `cond_*` | built-in | thermal conduction |
| `dhj_*`, `dhjsci_*`, `dhjph_*`, `dhjmhd_*`, `dhjnp_*`, `dhjrst_*`, `dhjmu1_*` | `deep_hot_jupiter_rt` | the science configs: spherical polar, etotgrav, user BCs, RT source, polar boundary, MHD, restart |
| `hb_*`, `hbt_*`, `hbg_*` | `hotbubble` | well-balanced source, limiter-sensitivity reproducer |
| `cv_*` | `convection` | etotgrav + `wellbalance_kappeli` |
| `tabsod_*` | built-in `shock_tube` | the TABULATED EOS end to end: cgs shock tube through the H ionization zone, run with `hydro/eos=general` (the input already sets `general_eos=table`). Since `4a59a1c9` its `pl`/`pr` are real PRESSURES (6.5e5 / 6.5e4), not the old energy proxies 1.144e6 / 1.144e5 |
| `sc_*` | `solar_convection` | (p,T) atmosphere; NEEDS the hardcoded `Rgas` patched to `83144621.4563013` to compare — see [[general-eos-project]] |

**Comparison scripts matter more than they look.** `.bin` output is SINGLE precision
(~6e-8 noise floor); `.tab` with `data_format = %24.16e` is double. Prefer `.tab` for
anything at round-off level. And NEVER normalise a difference by a per-variable max that
can be near zero — velocities against the sound speed, fields against max|B|, and always
look at the absolute difference too.

`eos_compare.py` is the driver that runs a case twice from ONE input with only
`<fluid>/eos` overridden — prefer it over hand-running pairs. It knows the per-case
overrides, including the `dhj` Rgas fix below.

**`dhj_gen.athinput` is the exception that still carries the physical `Rgas = 4.593e7`**, so
running it as-is shows an ORDER-UNITY ideal-vs-general difference that is purely the mu
artifact. Override `problem/Rgas=83144621.4563013`; it then agrees exactly.

Note `dhjsci_*` etc. use `Rgas = 83144621.4563013` (k_B/m_u, i.e. mu = 1) rather than the
physical `4.593e7`, because the gamma-law placeholder carries mu = 1 and the comparison is
otherwise off by that factor.
