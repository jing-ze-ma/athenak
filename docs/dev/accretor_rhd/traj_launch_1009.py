# central ballistic trajectory L1 -> r_out (and -> photosphere r 9.39) vs launch speed (pgen launch: x_L1 - 1e-4, v = -c/(a Omega) x)
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
G,MSUN,RSUN=6.674e-8,1.989e33,6.957e10; GU=G*MSUN/(RSUN*1e10)
MA,MD=16.,18.; P=3.81
a=(G*34*MSUN*(P*86400)**2/(4*np.pi**2))**(1/3)/RSUN; om=np.sqrt(GU*34/a**3); aom=a*om
mu=MD/34; xa,xd=-mu,1-mu
gx=lambda x: -(1-mu)*(x-xa)/abs(x-xa)**3-mu*(x-xd)/abs(x-xd)**3+x
xl1=brentq(gx,xa+1e-3,xd-1e-3)
def rhs(t,s):
    x,y,vx,vy=s; r1=np.hypot(x-xa,y); r2=np.hypot(x-xd,y)
    return [vx,vy,-(1-mu)*(x-xa)/r1**3-mu*(x-xd)/r2**3+x+2*vy,-(1-mu)*y/r1**3-mu*y/r2**3+y-2*vx]
for c in (18.6791,45.45,69.7):
    for rc in (13.7926,9.39246):
        ev=lambda t,s: np.hypot(s[0]-xa,s[1])-rc/a; ev.terminal=True; ev.direction=-1
        so=solve_ivp(rhs,[0,6],[xl1-1e-4,0,-c/aom,0],events=ev,rtol=1e-11,atol=1e-13)
        s=so.y_events[0][0]; dx,y=s[0]-xa,s[1]; r=np.hypot(dx,y)
        vr=(s[2]*dx+s[3]*y)/r*aom; vp=(-s[2]*y+s[3]*dx)/r*aom
        print('launch %.2f km/s r %.4f: phi %.3f deg v_r %.2f v_phi %.2f |v| %.2f t %.4f code'%(c,rc,np.degrees(np.arctan2(y,dx)),vr,vp,np.hypot(vr,vp),so.t_events[0][0]/om))
# transverse curvature at the r_out crossing of the launch-c_T(L1) trajectory (as plaskett_setup sec. 2)
GMA,GMD=GU*MA,GU*MD; XCM=a*MD/34
def roche3(x,y,z):
    return -GMA/np.sqrt(x*x+y*y+z*z)-GMD/np.sqrt((x-a)**2+y*y+z*z)-0.5*om**2*((x-XCM)**2+y*y)
for c in (18.6791,45.45):
    ev=lambda t,s: np.hypot(s[0]-xa,s[1])-13.7926/a; ev.terminal=True; ev.direction=-1
    so=solve_ivp(rhs,[0,6],[xl1-1e-4,0,-c/aom,0],events=ev,rtol=1e-11,atol=1e-13)
    s=so.y_events[0][0]; xo,yo=(s[0]-xa)*a,s[1]*a
    vh=np.array([s[2],s[3]]); vh/=np.hypot(*vh); nh=np.array([-vh[1],vh[0]]); h=1e-3
    f0=roche3(xo,yo,0)
    Pnn=(roche3(xo+h*nh[0],yo+h*nh[1],0)-2*f0+roche3(xo-h*nh[0],yo-h*nh[1],0))/h**2/om**2
    Pzz=(roche3(xo,yo,h)-2*f0+roche3(xo,yo,-h))/h**2/om**2
    print('launch %.2f: Phi_nn %.4f Phi_zz %.4f Omega^2 at the r_out crossing'%(c,Pnn,Pzz))
