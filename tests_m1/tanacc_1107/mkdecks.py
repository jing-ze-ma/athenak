"""tanacc-1107: decks for the full-tensor (vet_gd_tangential) accuracy test.
geometry g = sp | cs, resolution lev = lo | hi, case = star | shd | semi | hot, tan = on | off,
dirs = geo2 (production 42 x vet_rotate_avg 8) | geo6 (362 x 8)."""
import math, os, sys
D = '/viper/ptmp2/jinma/openitems_1107/tanacc/inp'
RES = {'sp': {'lo': (32, 36, 48), 'hi': (64, 72, 96)},
       'cs': {'lo': (32, 48, 48), 'hi': (64, 96, 96)}}
CASES = {  # clump_rho (chi = 1e-30 rho), clump_temp
    'star': None,
    'shd': (1.0e33, 1.0e-30),            # tau_diam 800, cold
    'semi': (1.25e30, 1.0e-30),          # tau_diam 1.0, cold absorber
    'mid': (2.0e31, 1.0e-30),            # tau_diam 16, cold absorber (VET faces)
    'hot5': (5.0e30, (0.1/1e16)**0.25),  # tau_diam 4, emitting, a T^4 = 0.1
    'hot': (5.0e31, (0.1/1e16)**0.25),   # tau_diam 40, a T^4 = 0.1 (arad 1e16: slow sound)
}
def deck(g, lev, case, tan, dirs='geo2'):
    n1, n2, n3 = RES[g][lev]
    if g == 'sp':
        mesh = f"""nx2 = {n2}
x2min = {math.pi/2 - 0.6!r}
x2max = {math.pi/2 + 0.6!r}
ix2_bc = periodic
ox2_bc = periodic
nx3 = {n3}
x3min = 0.0
x3max = {math.pi/2!r}
ix3_bc = periodic
ox3_bc = periodic
use_spherical_polar = true"""
        mb = (n1, n2, n3)
        cph = math.pi/4
        hx = "sp_x2_periodic_image = true\n"
        wall = "vet_gd_wall_interp = true\n"
    else:
        mesh = f"""nx2 = {n2}
x2min = -1.0
x2max = 1.0
ix2_bc = panel
ox2_bc = panel
nx3 = {n3}
x3min = -1.0
x3max = 1.0
ix3_bc = panel
ox3_bc = panel
use_cubed_sphere = true"""
        mb = (n1, n2, n3)
        cph = 0.0
        hx = ""
        wall = ""
    if dirs == 'geo2':
        dset = "vet_gd_gl_nmu = 0\nvet_gd_geo = 2"
    elif g == 'sp':   # high-direction arm on the wedge: GL 18 x 24 (wall-invariant)
        dset = "vet_gd_gl_nmu = 18\nvet_gd_gl_nphi = 24\nvet_gd_gl_stagger = true"
    else:             # high-direction arm on cs: geodesic level 6 (362)
        dset = "vet_gd_gl_nmu = 0\nvet_gd_geo = 6"
    cl = ""
    if CASES[case] is not None:
        rho, t = CASES[case]
        cl = f"""clump_r = 3.0
clump_th = {math.pi/2!r}
clump_ph = {cph!r}
clump_a = 0.4
clump_rho = {rho!r}
clump_temp = {t!r}
"""
    return f"""<comment>
problem = tanacc-1107 {g} {lev} {case} tan={tan} {dirs}: black-body sphere r=1 (F_in = 1) in a
          thin medium + clump at r = 3 (a = 0.4); exact = ray quadrature (ana/exact.py)

<job>
basename = ta

<mesh>
nghost = 2
nx1 = {n1}
x1min = 1.0
x1max = 6.0
ix1_bc = outflow
ox1_bc = outflow
{mesh}

<meshblock>
nx1 = {mb[0]}
nx2 = {mb[1]}
nx3 = {mb[2]}

<time>
evolution = dynamic
integrator = rk2
cfl_number = 0.3
nlim = 60
tlim = 1e+30
ndiag = 10

<hydro>
eos = ideal
reconstruct = plm
rsolver = hllc
gamma = 1.666666666666667
rad_signal_speed = false
{hx}
<rad_m1>
c_light = 1.0
chat_over_c = 1.0
e_floor = 1.0e-30
opacity = const
kappa_p = 1.0e-30
kappa_e = 1.0e-30
kappa_f = 1.0e-30
kappa_s = 0
arad = 1.0e16
coupling = true
gas_feedback = false
transport = implicit
time_scheme = be
implicit_cfl = 1000.0
implicit_solver = bicgstab
implicit_precond = mg
implicit_maxit = 200
implicit_tol = 1.0e-10
implicit_lin_tol = 1.0e-12
implicit_closure_lag = step
marshak_q = 0.5
ic_rad = none
implicit_bc_x1min = flux
implicit_flux_x1min = 1.0
implicit_bc_x1max = marshak
implicit_ebath_x1max = 0.0
vet_bc_inner = emit
closure = vet_col
vet_gd = true
vet_gd_replace = true
vet_col_lat_taucut = 30.0
{dset}
vet_rotate_avg = 8
vet_gd_rotate_every = 2
vet_col_lat_offdiag = operator
implicit_op_stencil = true
vet_gd_tangential = {'true' if tan == 'on' else 'false'}
{wall}vet_gd_twin = false
vet_gd_twin_fuse = false
vet_gd_symref = none
vet_source_noesrc = false
vet_scatter = true
vet_scatter_form = ma
vet_scatter_j = relaxed
implicit_face_opac_n = true
implicit_pos_floor = true
implicit_pos_floor_solve = true
vet_gd_allow_clamp = true
implicit_flux = blend
implicit_flux_faces = all
implicit_flux_beam = halfrange
implicit_blend = idort_f
implicit_blend_tau0 = 1
implicit_blend_alpha = 1
implicit_blend_flo = 0.3
implicit_blend_fhi = 0.6
implicit_blend_xthin = 0
implicit_thin_switch = fse
implicit_thin_corr = hc
implicit_thin_js_face = hcface
implicit_thin_hcface_gd = x1
implicit_thin_hcface = tau
implicit_thin_bq_gd = cell
implicit_thin_js_gd = x1
implicit_thin_hcsf = j
implicit_thin_recon = jse
implicit_thin_bq = lin
implicit_hr_face_gd = none
vet_interp = auto
vet_col_lat_cut_thin = 1.0
vet_seg = jump
vet_bc_diff = face

<problem>
pgen_name = rad_m1_beam
m1_test = sph_shell
gas_rho = 1.0e24
gas_temp = 1.0e-30
e_out = 1.0e-8
e_amp = 0.0
r0 = 2.0
width = 0.2
{cl}gas_hold = true
user_srcs = true

<output1>
file_type = bin
variable = m1
dt = 1e+30

<output2>
file_type = bin
variable = m1_lat
dt = 1e+30

<output3>
file_type = bin
variable = m1_face
dt = 1e+30
"""
def extra():
    names = []
    for g in ('sp', 'cs'):
        for lev in ('lo', 'hi'):
            for case in ('mid', 'hot5'):
                for tan in ('on', 'off'):
                    nm = f'{g}_{lev}_{case}_{tan}'
                    open(f'{D}/{nm}.athinput', 'w').write(deck(g, lev, case, tan))
                    names.append(nm)
    open(D + '/LIST2', 'w').write('\n'.join(names) + '\n')


if __name__ == '__main__' and len(sys.argv) > 1:
    extra()
elif __name__ == '__main__':
    os.makedirs(D, exist_ok=True)
    names = []
    for g in ('sp', 'cs'):
        for lev in ('lo', 'hi'):
            for case in CASES:
                for tan in ('on', 'off'):
                    nm = f'{g}_{lev}_{case}_{tan}'
                    open(f'{D}/{nm}.athinput', 'w').write(deck(g, lev, case, tan))
                    names.append(nm)
        for case in ('shd', 'hot'):
            for tan in ('on', 'off'):
                nm = f'{g}_lo_{case}_{tan}_hd'
                open(f'{D}/{nm}.athinput', 'w').write(deck(g, 'lo', case, tan, 'hd'))
                names.append(nm)
    print(len(names)); open(D + '/LIST', 'w').write('\n'.join(names) + '\n')
