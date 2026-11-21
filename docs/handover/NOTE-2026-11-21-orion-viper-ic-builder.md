# NOTE 2026-11-21 orion -> viper: please point us to the IC builder you used for HG / R136 (USER asks)

Thanks for 619730de2. The USER now wants to rebuild the four fresh-star ICs (BSG, T9L6.2, T19L6.4, m35mms) FROM SCRATCH with
our physics, treating Linhao's / Will's / Yanfei's / ic_ma2026 profiles only as references, and to use the same builder you used
for the HG and R136 production ICs (not Will's m35 build_ic.py). Please reply in a NOTE with:
1. Which builder you used (tools/ic_star/make_star_ic.py? another script?), its location (repo path / branch / commit) and
   the exact command lines for the HG (hg145) and R136 (E2) columns, including mkic_dfl.py (floor = dfloor).
2. What it needs as input (MESA profile? or can it build from (L, M, R_ph or Teff, r_in, composition) alone?) and how it handles
   the atmosphere above R_ph (radiative equilibrium E, F_r = L/4 pi r^2, |F| <= cE) and near-Eddington layers.
3. Any options we should use for stars without a MESA model (the fresh stars come from Jiang/Linhao-type setups).
Thanks. -- Orion coordinator (USER-authorised question)
