"""dynamo.py (wasp189b): dynamo field from the energy-flux scaling (Christensen+2009;
Reiners & Christensen 2010 eq. 1, B_dyn = 4.8e3 G (M L^2 / R^7)^(1/6), solar units; as used for
hot Jupiters by Yadav & Thorngren 2017), L = 4 pi R^2 sigma T_int^4 with the RCE T_int
(Thorngren+2019 fit as in rce_w121.py).  The surface polar dipole is scaled from Jupiter
(B_dip,pole 8.3 G from its dipole moment; Jupiter T_int 99 K) by the ratio of B_dyn."""
import numpy as np
SIG, RSUN, LSUN, MJ_MSUN, RJ = 5.670374419e-5, 6.957e10, 3.828e33, 9.5458e-4, 7.1492e9


def bdyn(m_mj, r_rj, tint):
    r = r_rj*RJ
    L = 4*np.pi*r**2*SIG*tint**4/LSUN
    return 4.8e3*(m_mj*MJ_MSUN*L**2/(r/RSUN)**7)**(1/6), L


bj, _ = bdyn(1.0, 1.0, 99.0)
for name, m, r, t in (('Jupiter', 1.0, 1.0, 99.0), ('WASP-121b (3x RCE)', 1.17, 1.742, 567.4),
                      ('WASP-189b (RCE T_int)', 1.99, 1.619, 490.5),
                      ('WASP-189b T_int 400 K', 1.99, 1.619, 400.0),
                      ('WASP-189b T_int 600 K', 1.99, 1.619, 600.0)):
    b, L = bdyn(m, r, t)
    print('%-24s M %.2f M_J R %.3f R_J T_int %5.1f K: L_int %.2e L_sun  B_dyn %6.1f G  '
          'B_dip,pole (Jupiter-scaled) %5.1f G' % (name, m, r, t, L, b, 8.3*b/bj))
