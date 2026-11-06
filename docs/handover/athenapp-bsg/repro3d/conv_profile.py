#!/usr/bin/env python3
# BSG/MESA_profile_combined_cgs.txt -> MESA_profile_combined.txt in his deck's code units
import numpy as np
ru,Tu,lu,mu=6.15245e-5,1.07887e6,6.957e10,.62
R=8.314462618e7/mu; a=7.565733250033928e-15; c=2.99792458e10
d=np.loadtxt('MESA_profile_combined_cgs.txt')
out=np.column_stack([d[:,0]/lu,d[:,1]/Tu,d[:,2]/ru,d[:,3]/(ru*lu**3),d[:,4]*ru*lu,d[:,5]/(c*a*Tu**4)])
np.savetxt('MESA_profile_combined.txt',out,fmt='%.12e')
i=np.argmin(abs(d[:,0]/lu-20)); print('r=20: T,rho code (expect 0.998, 0.993):',out[i,1],out[i,2])
