import numpy as np
import json
d = np.loadtxt('/viper/u2/jinma/ATHENAK/bench/m1_stage2/ic/ic_m1_V2.txt')
z, rho, T, Pg, E, F, f, chi, kR, kP, tau, geff = [d[:, i] for i in range(12)]
g0 = 3.98107e5
arad = g0 - geff
c = 2.99792458e10
a_cgs = 7.5657e-15
Prad = E / 3.0
beta = Pg / (Pg + Prad)
# gradients (centred)
lnT = np.log(T)
lnPg = np.log(Pg)
nab_gas = np.gradient(lnT, lnPg)
lnPt = np.log(Pg + Prad)
nab_tot = np.gradient(lnT, lnPt)
Hp_gas = Pg / (rho * geff)              # gas-pressure scale height under g_eff
Hp_tot = (Pg + Prad) / (rho * g0)
nad = 0.4
N2 = (geff / Hp_gas) * (nad - nab_gas)
sig = np.sqrt(np.maximum(-N2, 0.0))
# mixture ad gradient (Chandrasekhar) for the true criterion
G1 = beta + (4 - 3 * beta)**2 * (1.0 - beta) * (4 + beta) / \
    (beta + 12 * (1 - beta) * (4 + beta))  # approx not used
nad_mix = (2.0 + 2.0 * (4.0 + beta) * (1.0 - beta) / beta**2) / \
    (32.0 - 24.0 * beta - 3.0 * beta**2) * 0 + (1.0)  # placeholder
# standard: nabla_ad = (Gamma2-1)/Gamma2 ; use Chandrasekhar formulas
chi_t = (4.0 - 3.0 * beta)  # dlnP/dlnT at const rho for mixture (gas+rad), P= Pg+Prad
# Chandrasekhar: Gamma_1 = beta + (4-3beta)^2 (gamma-1)/(beta+12(gamma-1)(1-beta))
gam = 5.0 / 3.0
G1c = beta + (4 - 3 * beta)**2 * (gam - 1) / (beta + 12 * (gam - 1) * (1 - beta))
G3m1 = (G1c - beta) / (4 - 3 * beta)
nad_mix = G3m1 / G1c
N2_mix = (g0 / Hp_tot) * (nad_mix - nab_tot)
sig_mix = np.sqrt(np.maximum(-N2_mix, 0.0))
print(" z          tau     beta  g/geff  nab_gas nad_g  sig_gas"
      "  nab_tot nad_mix sig_mix  arad/g")
for i in range(0, len(z), 3):
    print(
        "%10.3e %8.3f %5.3f %6.2f %7.3f %5.2f %8.4f %7.3f %7.3f %8.4f %6.3f" %
        (z[i],
         tau[i],
         beta[i],
         g0 /
         geff[i],
         nab_gas[i],
         nad,
         sig[i],
         nab_tot[i],
         nad_mix[i],
         sig_mix[i],
         arad[i] /
         g0))
json.dump({'z': z.tolist(),
           'tau': tau.tolist(),
           'beta': beta.tolist(),
           'g_over_geff': (g0 / geff).tolist(),
           'nabla_gas': nab_gas.tolist(),
           'sigma_gas': sig.tolist(),
           'nabla_tot': nab_tot.tolist(),
           'nabla_ad_mix': nad_mix.tolist(),
           'sigma_mix': sig_mix.tolist(),
           'arad_over_g': (arad / g0).tolist(),
           'Hp_gas': Hp_gas.tolist()},
          open('h1.json',
               'w'))
print("max sigma_gas", sig.max(), "at z", z[np.argmax(sig)], "tau", tau[np.argmax(sig)])
print("max sigma_mix", sig_mix.max(), "at z",
      z[np.argmax(sig_mix)], "tau", tau[np.argmax(sig_mix)])
