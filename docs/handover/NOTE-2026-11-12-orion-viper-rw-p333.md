# NOTE 2026-11-12 Orion -> viper: rw-p333 (P 3.33 tau 10, why the order is lost)

FINAL binary 84cd20888 (Release, MPI, -ffp-contract=off; the rw-amp1e6 binary). One rank per case, host orion003.
Wall: cases_base 300 s, cases_diag 157 s, cases_a 306 s (37 cores). The pre-check passed: dt_retry false, implicit_tol 1e-12, tol_rel 0,
cnorm 1e-13, amp 1e-6, outputs m1_vet + m1_vsc.

**Harness note.** All 29 cases_a tags are 240+ characters with the long arm keys, and the first launch failed on Errno 36 (file name too long, GPFS 255 limit).
I patched rw_amp1e6/run_radwave.py on Orion: case_tag() hashes tags over 250 characters to `tag[:120] + "_h" + md5[:12]` (shorter tags unchanged, so
every viper dir name is still valid). The rerun used a 200-character threshold for its directory names, so the 8 ladder cases sit in hashed dirs with full-name symlinks. Nothing else changed.

**Summary.**
- (a) dt and (b) dx both converge at 2nd order to the same **omega_inf ~ (13.34565, -3.66302)**, which is off vetqs_d by d w_im = -2.39e-3.
- (d) Of the 7 arms, none moves omega_inf except physflux=false. That arm is a different model: 13.62 - 3.77i, flat in N, so not a discretisation error.
- implicit_closure_lag=pass and vet_scatter=false are bitwise equal to the baseline. vet_source=gas and sc=rebuild move omega by 1e-5. upwind enthalpy is 1st order to the same limit.
- implicit_flux=central (needs extra keys, see below) lands on the baseline to 1e-5.
- (c) The code's VET perturbation matches the vetqs_d closure: R = E0 dD/(qT T1) = 0.99981 (N128) -> 0.99995 (N256), i.e. -> 1 at 2nd order.
- **Decisive cross-check (rw-edd NOTE 25e404e72).** With closure = eddington and its exact analytic root, the same cell shows the same signature: omega_inf off by d w_im -2.8e-3, p 0.91 / 1.29.
  The other 17 cells are clean. So it is not a vetqs_d model mismatch. It is a term the code solves differently from the lab-frame linear system in
  radwave_disp, visible only at P 3.33 with tau ~ 10. qT x 0.99 / 1.01 (variants.py) moves w_im by -/+4.2e-3 but w_re the wrong way for the observed offset, so it is not a pure qT rescale.

## 1a. python3 reana.py cases_base.txt runs
```
  16    32 2.354 13.04719287 -3.69285511 1.1829e-02 1.5112e-02 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc
  32    64 2.354 13.26817051 -3.66615980 2.8865e-03 3.9549e-03 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc
  64   128 2.354 13.32598583 -3.66333183 7.1950e-04 1.0236e-03 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc
 128   256 2.354 13.34069340 -3.66303897 1.8482e-04 2.7682e-04 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc
 256   512 2.354 13.34439795 -3.66301404 8.5125e-05 1.0730e-04 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc
1024  2048 2.354 13.34556182 -3.66301577 8.3400e-05 7.5888e-05 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc
```

## 1b. python3 reana.py cases_a.txt runs
```
 512  1024 2.354 13.34532732 -3.66301425 8.1911e-05 7.9645e-05 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc
 256  1024 1.177 13.34404684 -3.66296495 1.0009e-04 1.0319e-04 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc
 256  2048 0.588 13.34395974 -3.66295325 1.0377e-04 1.0231e-04 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc
 256  4096 0.294 13.34393926 -3.66295060 1.0464e-04 1.0213e-04 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc
  32  4096 0.037 13.24305799 -3.66061461 3.5220e-03 3.4654e-03 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc
  64  4096 0.074 13.31902003 -3.66222496 8.8672e-04 9.2088e-04 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc
 128  4096 0.147 13.33886811 -3.66277917 2.4315e-04 2.5488e-04 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc
 512  4096 0.588 13.34521966 -3.66299990 8.6175e-05 7.8633e-05 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc
 128   256 MISSING expected non-empty vector for x
 256   512 MISSING expected non-empty vector for x
 512  1024 MISSING expected non-empty vector for x
 128   256 2.354 13.34055965 -3.66286603 1.8994e-04 2.6484e-04 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc rm_vet_source=gas
 256   512 2.354 13.34436395 -3.66297485 8.5702e-05 1.0405e-04 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc rm_vet_source=gas
 512  1024 2.354 13.34531863 -3.66300565 8.1946e-05 7.8902e-05 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc rm_vet_source=gas
 128   256 2.354 13.34070097 -3.66313288 1.8468e-04 2.7608e-04 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc rm_time2_vet_sc=rebuild
 256   512 2.354 13.34438562 -3.66302487 8.5938e-05 1.0799e-04 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc rm_time2_vet_sc=rebuild
 512  1024 2.354 13.34532197 -3.66301485 8.2025e-05 7.9828e-05 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc rm_time2_vet_sc=rebuild
 128   256 2.354 13.34069340 -3.66303897 1.8482e-04 2.7682e-04 t=0.470779 amp=1e-06 rm_implicit_closure_lag=pass rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc
 256   512 2.354 13.34439795 -3.66301404 8.5125e-05 1.0730e-04 t=0.470779 amp=1e-06 rm_implicit_closure_lag=pass rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc
 512  1024 2.354 13.34532732 -3.66301425 8.1911e-05 7.9645e-05 t=0.470779 amp=1e-06 rm_implicit_closure_lag=pass rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc
 128   256 2.354 13.34069340 -3.66303897 1.8482e-04 2.7682e-04 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc rm_vet_scatter=false
 256   512 2.354 13.34439795 -3.66301404 8.5125e-05 1.0730e-04 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc rm_vet_scatter=false
 512  1024 2.354 13.34532732 -3.66301425 8.1911e-05 7.9645e-05 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc rm_vet_scatter=false
 128   256 2.354 13.34305069 -3.66044282 2.8485e-03 2.9139e-03 t=0.470779 amp=1e-06 rm_implicit_enthalpy=upwind rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc
 256   512 2.354 13.34510412 -3.66227869 1.4062e-03 1.4521e-03 t=0.470779 amp=1e-06 rm_implicit_enthalpy=upwind rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc
 512  1024 2.354 13.34545063 -3.66283082 7.0243e-04 7.2790e-04 t=0.470779 amp=1e-06 rm_implicit_enthalpy=upwind rm_implicit_hr_cmem_rst=true rm_implicit_physflux_zold=fc
 128   256 2.354 13.61902822 -3.77311698 1.4101e-02 1.2966e-02 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux=false rm_implicit_physflux_zold=fc
 256   512 2.354 13.61990539 -3.77334528 1.4196e-02 1.3092e-02 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux=false rm_implicit_physflux_zold=fc
 512  1024 2.354 13.61952112 -3.77350940 1.4215e-02 1.3110e-02 t=0.470779 amp=1e-06 rm_implicit_hr_cmem_rst=true rm_implicit_physflux=false rm_implicit_physflux_zold=fc
```
The 3 MISSING rows are the arm implicit_flux=central (N 128/256/512): FATAL at init, `<rad_m1>/implicit_flux_faces = all needs implicit_flux = berthon | blend`.
Adding faces x1 + blend tau_f + beam closure FATALs in turn on `implicit_thin_switch needs implicit_flux = blend and implicit_flux_beam = halfrange`.
Adding thin_switch none FATALs on `implicit_physflux needs implicit_flux_beam = halfrange or an AP-HLL face`. Adding physflux=false runs:

python3 reana.py cases_central4.txt runs  (= central, faces x1, blend tau_f, beam closure, thin_switch none, physflux false; rw_p333/cases_central4.txt):
```
 128   256 2.354 13.34066888 -3.66302430 1.8569e-04 2.7744e-04 t=0.470779 amp=1e-06 rm_implicit_blend=tau_f rm_implicit_flux=central rm_implicit_flux_beam=closure rm_implicit_flux_faces=x1 rm_implicit_hr_cmem_rst=true rm_implicit_physflux=false rm_implicit_physflux_zold=fc rm_implicit_thin_switch=none
 256   512 2.354 13.34438422 -3.66300552 8.5291e-05 1.0746e-04 t=0.470779 amp=1e-06 rm_implicit_blend=tau_f rm_implicit_flux=central rm_implicit_flux_beam=closure rm_implicit_flux_faces=x1 rm_implicit_hr_cmem_rst=true rm_implicit_physflux=false rm_implicit_physflux_zold=fc rm_implicit_thin_switch=none
 512  1024 2.354 13.34532003 -3.66300969 8.1877e-05 7.9626e-05 t=0.470779 amp=1e-06 rm_implicit_blend=tau_f rm_implicit_flux=central rm_implicit_flux_beam=closure rm_implicit_flux_faces=x1 rm_implicit_hr_cmem_rst=true rm_implicit_physflux=false rm_implicit_physflux_zold=fc rm_implicit_thin_switch=none
```
Compare physflux=false on halfrange (cases_a rows above): 13.6195 - 3.7735i. So on central faces the physflux-free scheme gives the
baseline omega_inf. physflux is what makes the halfrange faces equal to central here; neither moves omega_inf off (13.3457, -3.6630).

Time ladder N256 (rows 2-4 + base N256 nt512): w_re 13.344398 / 13.344047 / 13.343960 / 13.343939 (bitwise the viper amp-1e-5 numbers).
Space at nt 4096: w_re 13.24306 / 13.31902 / 13.33887 / 13.34394 / 13.34522 (N 32..512); differences 7.6e-2, 2.0e-2, 5.1e-3, 1.3e-3 (ratio 3.8, 3.9, 4.0).
w_im -3.66061 / -3.66222 / -3.66278 / -3.66295 / -3.66300 (ratio 2.9, 3.2, 3.4); Richardson limit (13.34565, -3.66302). The base N1024 nt2048 is 13.34556 - 3.66302i.

## 2. python3 diag_ana.py runs/<N128 / N256 nt4096 diagvet>
```
=== python3 diag_ana.py vetf_impl_x1_P3.33333_t10_N128_nt4096_amp1e-06_armTIME_diagvet_refvetq...
t=0.02942  R=E0 dD/(qT T1)=0.999810-0.000001j  |  J1/Jref=0.999801-0.000001j  |  (E1/T1)code/eig=0.999997+0.000004j  (T1/rho1)code/eig=0.999910+0.000217j  D0=0.333333333
t=0.05885  R=E0 dD/(qT T1)=0.999810-0.000001j  |  J1/Jref=0.999801-0.000001j  |  (E1/T1)code/eig=0.999997+0.000004j  (T1/rho1)code/eig=0.999701+0.000393j  D0=0.333333333
t=0.08827  R=E0 dD/(qT T1)=0.999810-0.000001j  |  J1/Jref=0.999801-0.000001j  |  (E1/T1)code/eig=0.999996+0.000004j  (T1/rho1)code/eig=0.999428+0.000514j  D0=0.333333333
t=0.11769  R=E0 dD/(qT T1)=0.999810-0.000001j  |  J1/Jref=0.999801-0.000001j  |  (E1/T1)code/eig=0.999996+0.000004j  (T1/rho1)code/eig=0.999089+0.000548j  D0=0.333333333
t=0.14712  R=E0 dD/(qT T1)=0.999810-0.000001j  |  J1/Jref=0.999801-0.000001j  |  (E1/T1)code/eig=0.999997+0.000004j  (T1/rho1)code/eig=0.998712+0.000446j  D0=0.333333333
t=0.17654  R=E0 dD/(qT T1)=0.999810-0.000001j  |  J1/Jref=0.999801-0.000001j  |  (E1/T1)code/eig=0.999997+0.000004j  (T1/rho1)code/eig=0.998364+0.000170j  D0=0.333333333
t=0.20597  R=E0 dD/(qT T1)=0.999810-0.000001j  |  J1/Jref=0.999801-0.000001j  |  (E1/T1)code/eig=0.999997+0.000004j  (T1/rho1)code/eig=0.998140-0.000271j  D0=0.333333333
t=0.23539  R=E0 dD/(qT T1)=0.999810-0.000001j  |  J1/Jref=0.999801-0.000001j  |  (E1/T1)code/eig=0.999996+0.000005j  (T1/rho1)code/eig=0.998127-0.000810j  D0=0.333333333
t=0.26481  R=E0 dD/(qT T1)=0.999810-0.000001j  |  J1/Jref=0.999801-0.000001j  |  (E1/T1)code/eig=0.999996+0.000005j  (T1/rho1)code/eig=0.998367-0.001333j  D0=0.333333333
t=0.29424  R=E0 dD/(qT T1)=0.999810-0.000001j  |  J1/Jref=0.999801-0.000001j  |  (E1/T1)code/eig=0.999995+0.000004j  (T1/rho1)code/eig=0.998832-0.001714j  D0=0.333333333
t=0.32366  R=E0 dD/(qT T1)=0.999810-0.000001j  |  J1/Jref=0.999801-0.000001j  |  (E1/T1)code/eig=0.999994+0.000003j  (T1/rho1)code/eig=0.999430-0.001857j  D0=0.333333333
t=0.35308  R=E0 dD/(qT T1)=0.999810-0.000001j  |  J1/Jref=0.999801-0.000001j  |  (E1/T1)code/eig=0.999995+0.000002j  (T1/rho1)code/eig=1.000036-0.001730j  D0=0.333333333
t=0.38251  R=E0 dD/(qT T1)=0.999810-0.000001j  |  J1/Jref=0.999801-0.000001j  |  (E1/T1)code/eig=0.999995+0.000002j  (T1/rho1)code/eig=1.000529-0.001367j  D0=0.333333333
t=0.41193  R=E0 dD/(qT T1)=0.999810-0.000001j  |  J1/Jref=0.999801-0.000001j  |  (E1/T1)code/eig=0.999996+0.000001j  (T1/rho1)code/eig=1.000832-0.000848j  D0=0.333333333
t=0.44135  R=E0 dD/(qT T1)=0.999810-0.000001j  |  J1/Jref=0.999801-0.000001j  |  (E1/T1)code/eig=0.999997+0.000002j  (T1/rho1)code/eig=1.000916-0.000267j  D0=0.333333333
t=0.47078  R=E0 dD/(qT T1)=0.999809+0.000000j  |  J1/Jref=0.999800+0.000000j  |  (E1/T1)code/eig=0.999998+0.000003j  (T1/rho1)code/eig=1.000793+0.000298j  D0=0.333333333
=== python3 diag_ana.py vetf_impl_x1_P3.33333_t10_N256_nt4096_amp1e-06_armTIME_diagvet_refvetq...
t=0.02942  R=E0 dD/(qT T1)=0.999952-0.000001j  |  J1/Jref=0.999951-0.000001j  |  (E1/T1)code/eig=0.999998+0.000001j  (T1/rho1)code/eig=0.999929+0.000044j  D0=0.333333333
t=0.05885  R=E0 dD/(qT T1)=0.999952-0.000001j  |  J1/Jref=0.999951-0.000001j  |  (E1/T1)code/eig=0.999998+0.000001j  (T1/rho1)code/eig=0.999823+0.000057j  D0=0.333333333
t=0.08827  R=E0 dD/(qT T1)=0.999952-0.000001j  |  J1/Jref=0.999951-0.000001j  |  (E1/T1)code/eig=0.999998+0.000001j  (T1/rho1)code/eig=0.999713+0.000052j  D0=0.333333333
t=0.11769  R=E0 dD/(qT T1)=0.999952-0.000001j  |  J1/Jref=0.999951-0.000001j  |  (E1/T1)code/eig=0.999998+0.000001j  (T1/rho1)code/eig=0.999589+0.000021j  D0=0.333333333
t=0.14712  R=E0 dD/(qT T1)=0.999952-0.000001j  |  J1/Jref=0.999951-0.000001j  |  (E1/T1)code/eig=0.999998+0.000000j  (T1/rho1)code/eig=0.999457-0.000059j  D0=0.333333333
t=0.17654  R=E0 dD/(qT T1)=0.999952-0.000001j  |  J1/Jref=0.999951-0.000001j  |  (E1/T1)code/eig=0.999998+0.000000j  (T1/rho1)code/eig=0.999344-0.000209j  D0=0.333333333
t=0.20597  R=E0 dD/(qT T1)=0.999952-0.000001j  |  J1/Jref=0.999951-0.000001j  |  (E1/T1)code/eig=0.999998+0.000001j  (T1/rho1)code/eig=0.999296-0.000423j  D0=0.333333333
t=0.23539  R=E0 dD/(qT T1)=0.999952-0.000001j  |  J1/Jref=0.999951-0.000001j  |  (E1/T1)code/eig=0.999998+0.000001j  (T1/rho1)code/eig=0.999355-0.000666j  D0=0.333333333
t=0.26481  R=E0 dD/(qT T1)=0.999952-0.000001j  |  J1/Jref=0.999951-0.000001j  |  (E1/T1)code/eig=0.999997+0.000001j  (T1/rho1)code/eig=0.999534-0.000875j  D0=0.333333333
t=0.29424  R=E0 dD/(qT T1)=0.999952-0.000001j  |  J1/Jref=0.999951-0.000001j  |  (E1/T1)code/eig=0.999997+0.000001j  (T1/rho1)code/eig=0.999807-0.000985j  D0=0.333333333
t=0.32366  R=E0 dD/(qT T1)=0.999952-0.000001j  |  J1/Jref=0.999951-0.000001j  |  (E1/T1)code/eig=0.999997+0.000000j  (T1/rho1)code/eig=1.000112-0.000956j  D0=0.333333333
t=0.35308  R=E0 dD/(qT T1)=0.999952-0.000001j  |  J1/Jref=0.999951-0.000001j  |  (E1/T1)code/eig=0.999997-0.000000j  (T1/rho1)code/eig=1.000375-0.000789j  D0=0.333333333
t=0.38251  R=E0 dD/(qT T1)=0.999952-0.000001j  |  J1/Jref=0.999951-0.000001j  |  (E1/T1)code/eig=0.999998-0.000001j  (T1/rho1)code/eig=1.000534-0.000525j  D0=0.333333333
t=0.41193  R=E0 dD/(qT T1)=0.999952-0.000001j  |  J1/Jref=0.999951-0.000001j  |  (E1/T1)code/eig=0.999998-0.000001j  (T1/rho1)code/eig=1.000564-0.000227j  D0=0.333333333
t=0.44135  R=E0 dD/(qT T1)=0.999952-0.000001j  |  J1/Jref=0.999951-0.000001j  |  (E1/T1)code/eig=0.999999-0.000000j  (T1/rho1)code/eig=1.000476+0.000045j  D0=0.333333333
t=0.47078  R=E0 dD/(qT T1)=0.999950-0.000000j  |  J1/Jref=0.999950+0.000000j  |  (E1/T1)code/eig=0.999999+0.000000j  (T1/rho1)code/eig=1.000302+0.000253j  D0=0.333333333
```

## 3. python3 variants.py
```
vetqs_d              (13.34636646957754-3.6606312478391607j)
vetqs_d + Doppler    (13.12721031617148-4.193944219126031j)
qT x 0.99           (13.35000805712874-3.6648358967006534j)
qT x 1.01           (13.342722978160433-3.6564276684891643j)
eddington            (13.701319925124452-4.086365212039564j)
qT (-1.0459274699443273-1.734723475976807e-16j) qv (-1.1564823173178715e-16-0.16646452695723907j)
```

Run trees: Orion /orion/ptmp/jinma/rwp333_1112/rw_p333 (runs/, runs_stale/ = empty dirs moved aside).
