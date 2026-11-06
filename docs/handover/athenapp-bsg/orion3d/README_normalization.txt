Orion Athena++ prod_seeded: normalization of bsg.hst and all files here.
bsg.hst is a copy of the run's history with its header: pi/3 sr wedge, code units, no factor applied.
Units (deck + nr_radiation/radiation.cpp, r_ideal=8.314462618e7/mu):
 density 6.15245e-05 g/cc, T 1.07887e+06 K, length 6.95700e+10 cm, mu 0.62
 v_unit=sqrt(R T/mu)=1.202834e+07 cm/s, t_unit=5783.83902 s = 6.694258e-02 d
 mass [3]: x rho_unit l^3 = 2.07164e+28 g
 energy [7-10]: x P_unit l^3 = 2.99727e+42 erg (P_unit=8.90143e+09 erg/cc); [10] is gas E (internal+KE), no Er, no gravity
 Er, Er0 [11,15]: x a_r T^4 l^3 = 3.45139e+42 erg (a_r T^4=1.02501e+10 erg/cc)
 Fr [12,16] are volume integrals of the flux: x c a_r T^4 l^3 = 1.03470e+53 erg cm/s. Local Fr x 3.07290e+20 erg/cm^2/s, L=4 pi r^2 Fr erg/s, full-sphere equivalent.
 GM = 260.795 code = 2.62502e+27 cgs = 19.780 Msun
 prat, crat as in the code: 1.15151e+00, 2492.38367
Viper's dE tables are full sphere (x12); dEgrav, dU are wedge. summary.txt: all wedge.
