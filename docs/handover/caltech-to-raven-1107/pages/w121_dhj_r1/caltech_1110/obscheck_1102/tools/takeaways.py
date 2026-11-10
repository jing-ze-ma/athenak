import json, numpy as np
O='/resnick/groups/carnegie_poc/jingze/obscheck_1102/'
B1=np.load(O+'out/bands_w121.npz',allow_pickle=True); B2=np.load(O+'out/bands_w76.npz',allow_pickle=True)
M=json.load(open(O+'page/obs_metrics.json'))
s=lambda B,b: B['st_'+b].item()
a,b,c,d=s(B1,'G395H_NRS1'),s(B1,'G395H_NRS2'),s(B1,'IRAC2'),s(B2,'IRAC2')
fe=M['fe']
T={'_intro':'Hydro p1 runs at rot 298-300 (5 dumps) post-processed with a 3-D GPU ray tracer (tools/em3d.py): rays through the cubed-sphere grid at 40 orbital phases, petitRADTRANS 3.4 R=1000 correlated-k (H2O, CO, CO2, CH4, OH, HCN, TiO, VO, FeH, Na, K, SiO) premixed by random overlap and rebinned to R~200, FastChem 3 equilibrium chemistry with local condensation (1x solar Asplund 2009 = the run EOS), H- bound-free/free-free, H2-H2/H2-He CIA, no clouds, no scattering, no reflected light; PHOENIX stars (Delrez 2016 / Ehrenreich 2020 parameters). Fe I: the fe_hires tool, per exposure. Model band = spread over the 5 dumps.',
'_order':['obs_summary.png','obs_w121_phasecurves.png','obs_w121_spectra.png','obs_w76_phasecurves.png','obs_fe_timeresolved.png'],
'obs_summary.png':dict(take='Not consistent: the model hotspot is ~%.0f deg east (observed 0-6 deg), the nightside is %.0fx too bright at 3-5 um and the dayside %.0f-%.0f %% too faint, in both planets.'%(a['off_max_deg'],b['night']/630e-6,100*(1-b['day']/4924e-6),100*(1-a['day']/3924e-6)),alt='model vs observed day flux, night flux, offset, brightness temperature per band',cap='Diamonds: model (Fourier fit of the phase curve; day = phase 0.5, night = phase 0). Dots: observed, 1 sigma. Brightness temperatures use the observed Rp/Rs and the same PHOENIX star.'),
'obs_w121_phasecurves.png':dict(take='WASP-121b: amplitude too small and peak too early in every band; NIRSpec NRS2 day %.0f vs 4924 ppm, night %.0f vs 630 ppm, peak offset %+.0f vs +2.7 deg.'%(b['day']*1e6,b['night']*1e6,b['off_max_deg']),alt='WASP-121b phase curves model vs JWST, Spitzer, HST',cap='Blue: model with the dump-to-dump spread. NIRSpec points: Evans-Soma 2025 band means (Zenodo). Spitzer: Morello 2023 and Davenport 2025. WFC3: Mikal-Evans 2022 binned. SOSS: Splinter 2025.'),
'obs_w121_spectra.png':dict(take='WASP-121b: the dayside spectrum has the right shape (1.4 um water and 4.5 um CO in emission, i.e. an inversion) but sits ~15-25 %% low; the nightside is 2-5x too bright.',alt='WASP-121b dayside and nightside emission spectra',cap='Data binned for display. Dayside: Kahle 2026 (NIRISS, NIRSpec at eclipse), Mikal-Evans 2020 (WFC3), Morello 2023 (Spitzer). Nightside: Evans-Soma 2025 phase 0.925-0.95, Mikal-Evans 2022, Morello 2023.'),
'obs_w76_phasecurves.png':dict(take='WASP-76b: Spitzer 4.5 um eclipse %.0f vs 3729 ppm, night %.0f vs ~800 ppm, peak offset %+.0f vs +0.7 deg; same pattern as WASP-121b.'%(d['day']*1e6,d['night']*1e6,d['off_max_deg']),alt='WASP-76b Spitzer phase curves and dayside spectrum',cap='May & Stevenson 2021 eclipse depth and derived nightside; Dang 2025; HST WFC3 Fu 2021 and Edwards 2020 (different dilution corrections). No JWST phase curve of WASP-76b is published.'),
'obs_fe_timeresolved.png':dict(take='Fe I through transit: WASP-76b matches Ehrenreich 2020 (model %.1f / %.1f km/s first / second half vs %.1f / %.1f, rms %.1f km/s); WASP-121b is too blueshifted late in transit (%.1f vs %.1f km/s).'%(fe['w76']['mod_first'],fe['w76']['mod_last'],fe['w76']['obs_first'],fe['w76']['obs_last'],fe['w76']['rms'],fe['w121']['mod_last'],fe['w121']['obs_last']),alt='Fe I line centre and depth vs orbital phase, model vs ESPRESSO',cap='Model: fe_hires per-exposure CCF with the pRT 3000 K Fe I template, planet rest frame, Gaussian fit; solid = both limbs as observed, dashed = one limb. Depths: CCF amplitude (template-weighted) vs observed CCF depth/contrast, so absolute depths are only roughly comparable. Data digitised where marked.')}
import os
if os.path.exists(O+'page/obs_clouds.json'):
    C=json.load(open(O+'page/obs_clouds.json'))
    def verdict(key, bn, nm):
        X=C[key][bn]; ob=X['obs']; cl=X['clear']; cs=X['cases']
        bn_=min(cs,key=lambda c: abs(c['night']-ob['night'][0])); bo=min(cs,key=lambda c: abs(c['off']-ob['off'][0])); bb=min(cs,key=lambda c:c['score'])
        return ('%s (%s): clouds can bring the nightside to %.0f ppm (observed %.0f, clear %.0f) and the offset to %+.0f deg (observed %+.1f, clear %+.0f); '
                'best overall case %s: night %.0f, amplitude %.0f vs %.0f, offset %+.0f deg, but the dayside stays at %.0f vs %.0f ppm (clouds do not form there).'%(
                nm, bn, bn_['night'], ob['night'][0], cl['night'], bo['off'], ob['off'][0], cl['off'], bb['name'], bb['night'], bb['amp'], ob['amp'][0], bb['off'], bb['day'], ob['day'][0]))
    cap='Post-processing only: condensate clouds placed in the frozen hydro T(p) (no radiative feedback, so the cloud does not heat the layers below or change the winds). Recipe (Parmentier+16/21, Roman & Rauscher 17/19, Roman+21, Frazier+26; notes in obscheck_1102/lit/clouds_notes.md): species condense where T < T_cond(p) (Visscher+10, Wakeford+17, 1x solar), cloud from the column base up 1 H, 3 H (0.43, 1.3 dex) or to the model top; species sets A = Al2O3+CaTiO3, C = A + Mg2SiO4+MgSiO3, B = C + Fe; full condensable mass (one f = 0.1 case); monodisperse radius 0.1/1/5 um, pRT Mie opacities; filled = extinction with (1-g)-scaled scattering counted as absorption, open = absorption only (no scattering in the tracer). Points: mean and spread over the 5 dumps, 20 phases each.'
    if 'w121' in C:
        T['obs_clouds_w121.png']=dict(take=verdict('w121','G395H_NRS2','WASP-121b'),alt='WASP-121b cloud grid: night flux, amplitude, offset, day flux vs particle radius',cap=cap)
    if 'w76' in C:
        T['obs_clouds_w76.png']=dict(take=verdict('w76','IRAC2','WASP-76b'),alt='WASP-76b cloud grid',cap=cap)
    b=C.get('best',{})
    parts=[]
    for k,v in b.items():
        ob=C[k][list(C[k].keys())[0]]['obs']
        okn=abs(v['night']-ob['night'][0])<max(0.3*ob['night'][0],3*ob['night'][1]); oko=abs(v['off']-ob['off'][0])<10; oka=abs(v['amp']-ob['amp'][0])<0.15*ob['amp'][0]
        parts.append('%s: night %s, amplitude %s, offset %s; dayside %.0f %% low'%({'w121':'WASP-121b','w76':'WASP-76b'}[k],'fixed' if okn else 'not fixed (%.0f vs %.0f ppm)'%(v['night'],ob['night'][0]),'fixed' if oka else 'not fixed (%.0f vs %.0f ppm)'%(v['amp'],ob['amp'][0]),'fixed' if oko else 'not fixed (%+.0f vs %+.1f deg)'%(v['off'],ob['off'][0]),100*(1-v['day']/ob['day'][0])))
    T['obs_clouds_summary.png']=dict(take='Best cloud deck per planet: '+'; '.join(parts)+'. Clouds do not form on the hot dayside, so they cannot remove the dayside deficit.',
        alt='best cloud-deck phase curve vs clear model vs data',cap='Orange: the cloud case with the smallest combined mismatch in night flux, amplitude and offset. '+cap)
    T['_order']=T['_order'][:4]+['obs_clouds_summary.png','obs_clouds_w121.png','obs_clouds_w76.png']+T['_order'][4:]
json.dump(T,open(O+'page/obs_takeaways.json','w'),indent=1)
print('\n'.join(v['take'] for k,v in T.items() if isinstance(v,dict)))
