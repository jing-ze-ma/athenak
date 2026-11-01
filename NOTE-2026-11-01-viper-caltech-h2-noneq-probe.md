# NOTE 2026-11-01 viper -> Caltech (USER): is equilibrium H/H2 valid on the WASP-121b night side? Yes below 1e-5 bar; frozen atomic above 1e-6 bar

FYI, nothing to run. A follow-up to NOTE-2026-11-01-viper-caltech-night-eta-floor-rain (651f6b733). After the x_e result, the user asked whether H/H2 also needs a relaxed tracer. We probed first.

## Method
- Script: viper /viper/ptmp2/jinma/h2probe_1101/probe.py. Data: the same s10 (10 G) rot-318 isobar stack as the night x_e estimate.
- Recombination time tau for local H atoms -> H2, starting all-atomic, as gas arriving from the day side would. Two channels:
  - three-body H + H + M: k3 = 6e-32 T^-0.25 + 2e-31 T^-0.5, M_eff = n_H + (n_H2 + n_He)/8;
  - H- + H -> H2 + e: k = 1.3e-9, with kinetic H- and the advected x_e (300x Saha).
- Flow time t_flow = pi r / |v_h| (1-2 days). Da = t_flow/tau.

## Result (night, zenith > 120 deg, median Da; the terminator is similar)

| log p [bar] | -8..-7 | -7..-6 | -6..-5 | -5..-4 | -4..-3 | -3..-2 | deeper |
|---|---|---|---|---|---|---|---|
| Da | 0.001 | 0.05 | 2 | 200 | 2e4 | 1e6 | > 1e8 |

- **p > 1e-5 bar: equilibrium H/H2 is valid.**
  - This covers the whole layer that carries heat by H2 dissociation and recombination.
  - It also covers the Fe-forming level, ~1e-4 bar (Da 60-200).
- **1e-6..1e-5 bar: marginal** (Da ~2).
- **p < 1e-6 bar: frozen atomic.**
  - Dayside H crosses the night side without recombining. The code instead puts equilibrium H2 there, with fa(eq) 0.01-0.07 at the night top.
  - So the night top is too heavy (mu ~2.3 vs ~1.3) and too compact, by up to ~1.7x in scale height.
  - That is the layer we already treat as a crude top.
- k3 x10 either way moves the boundary by ~0.5 dex in p.

## Decision (user): no H2 tracer; caveat only
- Unlike x_e, which enters only eta, f_H2 is part of the EOS (mu, internal energy incl. the 4.5 eV binding energy, gamma) and of the ck/CE tables, which all assume equilibrium. A non-equilibrium H2 would need a three-variable EOS and new tables.
- Possibly relevant to you: the night-top density and scale height bear on your Roche-lobe / dayside r(1e-6 bar) / one-way top questions.
  - The direction: real gas is puffier at the night top than the model.
  - If that layer ever matters, the cheap test is a run with the top composition pinned to atomic H.
