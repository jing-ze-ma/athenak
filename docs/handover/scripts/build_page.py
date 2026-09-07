import base64, numpy as np, json
S='/tmp/claude-28895/-viper-u2-jinma-ATHENAK-athenak/967431f9-6fb9-4e65-a716-2a3c3456d615/scratchpad'
def img(fn): return 'data:image/png;base64,'+base64.b64encode(open(S+'/'+fn,'rb').read()).decode()
out=np.load(S+'/jet_series.npy',allow_pickle=True).item()
def at(k,rot,col):
    a=out[k]; i=np.argmin(np.abs(a[:,0]-rot)); return a[i,col]
rows=[]
for rot in [8,12,16,20]:
    rows.append(f"<tr><td>{rot}</td><td>{at('cs hydro',rot,2):.2f}</td><td>{at('cs MHD (fixed BC, rot 8-20)',rot,2):.2f}</td><td>{at('cs MHD (drifting BC, rot 0-26)',rot,2):.2f}</td><td>{at('sp hydro',rot,2):.2f}</td></tr>")
tbl='\n'.join(rows)
hlast=out['cs hydro'][-1]; slast=out['sp hydro'][-1]
html=f"""<title>Deep Hot Jupiter on the Cube</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root{{--bg:#f3f5f7;--panel:#ffffff;--ink:#16202b;--muted:#5d6b7a;--line:#d6dde4;--hyd:#1f5fbf;--mhd:#c8451e;--sp:#5c5c5c;--accent:#0f6b5c}}
@media (prefers-color-scheme: dark){{:root:not([data-theme="light"]){{--bg:#12171d;--panel:#1a2129;--ink:#e7ecf1;--muted:#9aa7b4;--line:#2c3640;--hyd:#6fa3ec;--mhd:#f0855f;--sp:#b0b0b0;--accent:#5fc9b4}}}}
:root[data-theme="dark"]{{--bg:#12171d;--panel:#1a2129;--ink:#e7ecf1;--muted:#9aa7b4;--line:#2c3640;--hyd:#6fa3ec;--mhd:#f0855f;--sp:#b0b0b0;--accent:#5fc9b4}}
body{{background:var(--bg);color:var(--ink);font-family:'IBM Plex Sans',system-ui,sans-serif;font-size:16px;line-height:1.55;margin:0}}
main{{max-width:1180px;margin:0 auto;padding:40px 24px 80px}}
header{{border-bottom:1px solid var(--line);padding-bottom:20px;margin-bottom:28px}}
h1{{font-size:30px;font-weight:600;margin:0 0 6px;text-wrap:balance;letter-spacing:-0.01em}}
h2{{font-size:20px;font-weight:600;margin:40px 0 10px;text-wrap:balance}}
.eyebrow{{font-family:'IBM Plex Mono',monospace;font-size:12px;letter-spacing:.08em;text-transform:uppercase;color:var(--muted)}}
p{{max-width:72ch}}
.verdict{{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:14px;margin:22px 0 10px}}
.verdict div{{background:var(--panel);border:1px solid var(--line);border-radius:6px;padding:14px 16px}}
.verdict b{{display:block;font-size:13px;color:var(--muted);font-weight:500;text-transform:uppercase;letter-spacing:.06em;margin-bottom:6px}}
.verdict span{{font-size:17px;font-weight:600}}
figure{{margin:18px 0 8px;background:var(--panel);border:1px solid var(--line);border-radius:6px;padding:10px;overflow-x:auto}}
figure img{{max-width:100%;display:block}}
figcaption{{font-size:14px;color:var(--muted);padding:8px 4px 2px;max-width:none}}
table{{border-collapse:collapse;font-family:'IBM Plex Mono',monospace;font-size:14px;font-variant-numeric:tabular-nums;margin:12px 0}}
th,td{{border-bottom:1px solid var(--line);padding:6px 14px;text-align:right}} th{{color:var(--muted);font-weight:500}} td:first-child,th:first-child{{text-align:left}}
.k{{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:6px;vertical-align:middle}}
code{{font-family:'IBM Plex Mono',monospace;font-size:14px}}
ul{{max-width:72ch}}
</style>
<main>
<header>
<div class="eyebrow">AthenaK · deep_hot_jupiter_rt · 2026-09-06</div>
<h1>Deep Hot Jupiter on the Cube</h1>
<p>Does the cubed-sphere MHD run reproduce the 3D temperature structure and the equatorial jet, and does the field slow the jet? Three runs at rotation 20 on the same physics: <span class="k" style="background:var(--hyd)"></span>cs hydro (cs_prod_hyd_rot), <span class="k" style="background:var(--mhd)"></span>cs MHD with the fixed bottom boundary (cs_inflow_bcfix, 248f1b77), and <span class="k" style="background:var(--sp)"></span>spherical-polar hydro (sp_dhj_hyd) as the reference grid. Temperature and pressure come from the run's own tabulated EOS (inverted, checked to 5e-4 against the code's column dump); winds are projected from the gnomonic tangent basis onto zonal and meridional directions; the field is taken in the orthonormal frame the code stores.</p>
</header>

<div class="verdict">
<div><b>3D temperature</b><span>Matches spherical polar at every level</span></div>
<div><b>Equatorial jet</b><span>Present on all grids, but weak: ~1 km/s at 1e-3 to 1e-2 bar</span></div>
<div><b>Effect of B through rot 20</b><span>Not slowed. MHD jet sits inside the hydro run's own scatter; plasma beta &ge; 100 everywhere</span></div>
</div>

<h2>1. Temperature on isobars</h2>
<p>Day-night contrast, the substellar hot region, its eastward shift with depth and the small nightside equatorial hot spot at 0.1 bar (longitude ~140&deg;) appear identically on the cubed sphere and on spherical polar. Deep layers (1 bar) are horizontally uniform to ~100 K on all three. There is no imprint of the six cube panels or of the cube vertices in any map.</p>
<figure><img src="{img('fig1_Tmaps.png')}" alt="Temperature maps at four pressure levels for three runs"><figcaption>Fig. 1. Temperature on the 1e-3, 1e-2, 0.1 and 1 bar isobars, 5&deg; bins, longitude measured from the substellar point (positive = eastward, the direction of rotation). Columns: cs hydro, cs MHD (fixed boundary), sp hydro, all at rotation 20.</figcaption></figure>

<h2>2. Vertical structure and the zonal wind</h2>
<p>The T(p) profiles at the substellar, antistellar and both terminator points coincide between cs hydro and cs MHD to the line width and agree with spherical polar; all three converge to one deep adiabat below 0.1 bar. The zonal-mean zonal wind has the same structure on every grid: a narrow eastward equatorial jet between 1e-4 and 1e-2 bar, broad westward flow at mid latitudes above 1e-3 bar, and a quiescent interior below 0.1 bar (under 0.3 km/s everywhere). The day-to-night flow itself is strong, up to 15 km/s at the top, but it is divergent and cancels in the zonal mean: this ultra-hot, slowly rotating planet (3.5-day period) has only a weak superrotating jet.</p>
<figure><img src="{img('fig2_Tp_winds.png')}" alt="T(p) profiles and zonal-mean zonal wind"><figcaption>Fig. 2. Top left: T(p) at four longitudes. Contours: zonal-mean zonal wind (lat, p) for the three runs, red eastward. Bottom middle: equatorial (|lat| &lt; 12.5&deg;) zonal-mean wind versus pressure. Bottom right: maximum |u| on each isobar; cs hydro and cs MHD coincide, spherical polar is ~25 % higher at 0.01 to 0.1 bar, which is its 3x finer equatorial cell (0.84&deg; against 2.8&deg;), not a physics difference.</figcaption></figure>

<h2>3. Is the jet slowed by the field?</h2>
<p>Not through rotation 20. In the rotation-20 snapshot the MHD jet peaks higher than the hydro one (about 1.4 versus 0.9 km/s near 4e-3 bar), but the time series shows that the jet on every grid oscillates between roughly 0.5 and 1.1 km/s at 0.01 bar from rotation to rotation, and the MHD run stays inside that scatter; below 0.1 bar the runs are identical. The reason is in the field itself: plasma beta is at least 100 on the 0.01 bar isobar and above 1e4 deeper, so the Lorentz force is dynamically negligible at this stage. The magnetic field is being wound up by the flow (|B| follows the wind pattern at 0.01 bar, and the deep field is a stretched, mostly toroidal 10 to 300 G structure), but it has not yet fed back on the dynamics. The comparison below tracks the equatorial zonal-mean wind at 0.01 bar through time; the "drifting BC" series is the earlier cs MHD run whose bottom boundary leaked energy (fixed in 248f1b77) and is shown only as context.</p>
<table><tr><th>rotation</th><th>cs hydro</th><th>cs MHD fixed BC</th><th>cs MHD drifting BC</th><th>sp hydro</th></tr>{tbl}<tr><th colspan=5 style="text-align:left">equatorial zonal-mean u at 0.01 bar, km/s</th></tr></table>
<figure><img src="{img('fig3_jet_time.png')}" alt="Equatorial jet speed versus time"><figcaption>Fig. 3. Equatorial zonal-mean zonal wind (top) and maximum |u| (bottom) on four isobars versus time. cs hydro runs to rotation {hlast[0]:.0f} in this window, sp hydro to {slast[0]:.0f}; the fixed-boundary cs MHD run covers rotations 8 to 20.</figcaption></figure>

<h2>4. The field at rotation 20</h2>
<figure><img src="{img('fig4_field.png')}" alt="Magnetic field maps"><figcaption>Fig. 4. cs MHD (fixed boundary): |B|, radial field and plasma beta on the 0.01, 1, 10 and 100 bar isobars. The field is weakest under the substellar point at 0.01 bar where the flow diverges, and beta never drops below ~100.</figcaption></figure>

<h2>What to watch as cs_mhd_prod runs</h2>
<ul>
<li>The deep field grows on both grids (300 G at the bottom by rotation 20); the jet can only be slowed once beta approaches unity in the jet region, which has not happened by rotation 20. Repeat Fig. 2 and 3 every 20 rotations of the new production run.</li>
<li>The maximum wind speed at 0.01 to 0.1 bar is resolution-limited on the cubed sphere (Fig. 3, bottom): the sp reference is 25 % faster with a 3x finer equatorial cell. The zonal mean is not affected.</li>
<li>The nightside 0.1 bar equatorial hot spot at longitude ~140&deg; is common to all three runs and worth a physical explanation (converging day-night flow) rather than a numerical one.</li>
<li>The spherical-polar hydro reference is the old binary; the new sp MHD production died at rotation 10 from a polar-row instability that the cubed sphere does not have, so the cs run is currently the only usable MHD production.</li>
</ul>
<p class="eyebrow">Scripts: scratchpad/dhjcs.py (EOS inversion, cubed-sphere geometry, regridding), figs.py, fig3.py. Data: bench/cs_prod_hyd_rot, bench/cs_inflow_bcfix, bench/cs_prod_mhd_rot, bench/sp_dhj_hyd.</p>
</main>
"""
open(S+'/dhj_cube_report.html','w').write(html); print('written', len(html)//1024, 'kB')
