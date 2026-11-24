# The clean staggered implicit M1–VET radiation scheme in AthenaK

This note describes the grey radiation–hydrodynamics scheme selected by `<rad_m1>/scheme = clean`
in this AthenaK branch. It covers what is solved, how it is discretised in space and time, how
the implicit system is solved, which conservation and robustness properties hold, and where each
piece lives in the source. The note is written for users who need to run, check and cite the
scheme. Every statement was checked against the source of commit `512cf768a`. Where the code
does not settle a point, the text says so (see also §12).

Notation. $E$ is the radiation energy density and $\mathbf F$ the radiation energy flux, both in
the lab frame. $\rho$, $\mathbf v$, $e$ and $T$ are the gas density, velocity, internal energy
density and temperature, and $\boldsymbol\beta=\mathbf v/c$. $a$ is the radiation constant
(`<rad_m1>/arad`, in code units). Opacities are per unit mass: $\kappa_P$ (Planck, emission),
$\kappa_E$ (energy mean, absorption), $\kappa_F$ (flux mean) and $\kappa_s$ (scattering).
Products with the density are written $k_P=\rho\kappa_P$, $k_E=\rho\kappa_E$ and
$\chi=\rho(\kappa_F+\kappa_s)$ (the total extinction). $h$ is the length of an implicit stage:
$h=\gamma\,\Delta t$ in a two-stage step, $h=\Delta t$ in a backward-Euler step.

---

## 1. Overview

- **Moments.** Grey, two-moment ("M1-type") equations for $E$ and $\mathbf F$ in the
  mixed frame: lab-frame moments with comoving-frame opacities. The $O(v/c)$ terms are kept,
  and so are the $O(\beta^2)$ terms that are leading order in dynamic diffusion.
- **Closure.** A variable Eddington tensor (VET). The pressure tensor is
  $\mathsf P=\mathsf D\,E$, where $\mathsf D=\mathsf K/J$ comes from a formal solution of the
  grey transfer equation. The formal solution uses short characteristics along a fixed set of
  directions. The same rays also give the half-range fluxes and a non-locality measure, which
  shape the face fluxes.
- **Layout.** "Staggered": $E$ lives at cell centres and the normal component of $\mathbf F$
  lives on cell faces (one scalar per face).
- **Time.** IMEX Runge–Kutta. The hydrodynamics is explicit (AthenaK's Heun/`rk2`
  integrator, unchanged). Radiation transport and gas–radiation exchange are implicit, with an
  L-stable two-stage SDIRK tableau. Backward Euler is the fallback.
- **Solver.** Each implicit stage reduces to a scalar system for $E$: the face momentum
  equation is eliminated face by face. The system is solved by Newton/Picard passes. Each pass
  runs a right-preconditioned BiCGStab with a line-relaxation multigrid preconditioner and a
  per-cell Newton solve for the gas temperature. The face flux $\mathbf F$ is then recovered
  from its own line-tridiagonal solve.
- **Properties.** Total (gas + radiation) energy is conserved to round-off, apart from booked
  floors. The radiative momentum exchange is in flux form. $E\ge0$ holds by construction
  through M-matrix operators and Zalesak flux limiting. A restart is bitwise identical to a
  straight run. Global sums can be made independent of the MeshBlock layout and of the rank
  count.

The scheme needs `<rad_m1>/transport = implicit`; `scheme = clean` is the default with it.
It is the only implicit M1 solver in this build: other `scheme` values are FATAL.

---

## 2. Equations solved

### 2.1 Radiation moments (continuum form)

The module integrates the mixed-frame system below, written here with $\hat c=c$ (see §2.4):

$$
\frac{\partial E}{\partial t} + \nabla\!\cdot\mathbf F = -c\,G^0,\qquad
\frac{\partial \mathbf F}{\partial t} + c^2\,\nabla\!\cdot(\mathsf D E) = -c^2\,\mathbf G ,
$$

with the four-force evaluated in the comoving frame and transformed back:

$$
E_0 = (1+\beta^2)E - 2\boldsymbol\beta\cdot\mathbf F/c + \boldsymbol\beta\cdot\mathsf P\cdot\boldsymbol\beta,\qquad
\frac{\mathbf F_0}{c} = \frac{\mathbf F}{c} - \boldsymbol\beta E - \boldsymbol\beta\cdot\mathsf P,
$$
$$
g^0 = k_E E_0 - k_P\,a T^4,\qquad \mathbf g = \chi\,\mathbf F_0/c,\qquad
G^0 = g^0 + \boldsymbol\beta\cdot\mathbf g,\qquad \mathbf G = \mathbf g + \boldsymbol\beta\,g^0 .
$$

In the flux equation the relaxation term reads $-c\chi(\mathbf F-\mathbf a E)$ with the
**advective velocity** $\mathbf a=\mathbf v+\mathsf D\mathbf v$. This is the lab-frame flux of
radiation enthalpy carried by the gas. In the diffusion limit the steady flux is
$\mathbf F\simeq\mathbf a E-(c/\chi)\nabla\!\cdot(\mathsf D E)$.

The β-dependent parts of the source follow Ma, Pakmor, Justham & de Mink (2025), eqs. (56)–(57).
The emission carries the factor $1-\beta^2/2$. The energy equation gets the additional source
$h c\,k_E\,[(\tfrac12\beta^2-\boldsymbol\beta\cdot\mathbf a/c)\,E+2\boldsymbol\beta\cdot\mathbf F/c]$
(per stage), and the face momentum equation gets
$\theta\,c\,\boldsymbol\beta_f\,(h c\,k_P aT^4-h c\,k_E E)_f$ (see §4). In the implicit solve these
terms are evaluated at the current Newton iterate (`CleanRightHandSide`, `CleanFaceWeights`,
`CleanStagBetaE` in `rad_m1_clean.hpp`). The $(\mathbf v\cdot\mathbf F)\mathbf v/c^3$
and $\tfrac12\beta^2\kappa\mathbf F/c$ momentum corrections are not kept.

### 2.2 Coupling to the gas

The gas takes exactly what the radiation gives up:

$$
\frac{\partial(\rho e_{\rm tot})}{\partial t}\Big|_{\rm rad} = c\,G^0,\qquad
\frac{\partial(\rho\mathbf v)}{\partial t}\Big|_{\rm rad} = \mathbf G .
$$

The EOS is gas-only: radiation pressure enters the gas only through $\mathbf G$. Any AthenaK EOS
may be used, including tabulated (general) EOS. The temperature root takes $c_v$ from the EOS
(`CleanExchangeTemperature`).

In the discrete scheme the momentum the gas receives at each face is minus the change of the
radiation momentum across that face over the stage:

$$
\Delta(\rho\mathbf v)\cdot\hat n_f\;=\;-\Big[\frac{F_f-G_f}{c^2}+h\,\partial_n(\mathsf D E)_f\Big],
$$

shared between the two half cells in proportion to their extinction, $\chi_c/(\chi_L+\chi_R)$
(`CleanMomentumChange`). The energy the gas receives is the local exchange together with the
work of that force (§6.3, §7.1).

### 2.3 Opacities

`<rad_m1>/opacity = const` uses `kappa_p`, `kappa_e` (default `kappa_p`), `kappa_f` (default
`kappa_p`) and `kappa_s` (default 0). `opacity = table` reads stellar Rosseland and Planck tables
(`M1TableOpacities`, `rad_m1_opacity.hpp`), with
$\kappa_F=\kappa_R$ (scattering included), $\kappa_P=\kappa_E=$ the Planck mean and
$\kappa_s=0$. With `vet_scatter = true`, $\kappa_R$ is split into absorption
$\kappa_R-\kappa_e$ and electron scattering $\kappa_e$ (Ma et al. 2026, eq. 7). The split
changes only the source function of the rays (§3.1), not the moment equations. Within a stage
the opacities are re-evaluated at the Newton iterate's temperature during the first
`clean_stag_opac_passes` passes (default 3) and then frozen (`clean_stag_opac_update`,
`CleanStageOpacity`).

### 2.4 Speed of light

`<rad_m1>/chat_over_c` exists in the module, but **the clean scheme requires
$\hat c=c$**: `CleanInit` stops with a FATAL error otherwise. The reduced-speed-of-light
validity check of the module (`rsla_warn`, `rsla_force`) therefore does not apply to clean runs.

### 2.5 Units

`c_light`, `arad`, `temp_unit_kelvin`, `rho_unit_cgs` and `kappa_unit` are typed in **code
units**. At start-up the module prints each one next to the value that `<units>` implies and warns
when they disagree.

---

## 3. Closure: the variable Eddington tensor from rays

### 3.1 The formal solution

At each closure build the grey transfer equation

$$
\hat n\cdot\nabla I = \chi\,(S-I)
$$

is solved along a fixed set of straight-line directions $\hat n_d$ with weights $w_d$
($\sum w_d=1$). The source function is the current state's own, written in energy-density units:

$$
S = \epsilon\,aT^4 + (1-\epsilon)\,E,\qquad \epsilon=\min(\kappa_P\rho/\chi,\,1),
$$

or, with `vet_scatter`, the absorption/scattering mix of Ma et al. with the absorption share
capped at 1 (`clean_source_convex`, default on). The relaxed source forms (`vet_source = relaxed*`)
are **not** used by the clean scheme unless `clean_ray_source_relaxed = true`. This keeps the
closure independent of $\Delta t$ (`CleanRays`).

The sweep depends on the mesh:

| mesh | sweep | file / function |
|---|---|---|
| spherical polar (3-D wedge, 2-D axisymmetric, full sphere) | short characteristics shell by shell: inward rays top-down, then outward rays bottom-up. The upwind point is where the backward ray meets the neighbouring shell, with $\ln\chi$ and $\ln S$ bilinear in $(\theta,\phi)$ and the source linear in optical depth along the segment. | `rad_m1_vetgd.cpp` `VetGdBuild` |
| cubed sphere | the same shell sweep; the lateral reads cross the panel seams through an unfolded band (bilinear inside a panel, linear across a seam, barycentric at a cube vertex) | `rad_m1_vetgd_cs.cpp` |
| Cartesian | short-characteristics sweep (`vet_sc`); optionally the same direction sets (`vet_sc_set = gd`) | `rad_m1_vet.cpp` `VetShortChar` |

Directions are fixed in the Cartesian frame, so curvature enters only through the ray
positions and there is no angular redistribution. At the bottom the rays start from the
diffusion-limit intensity. At the top they start from vacuum (Marshak boundary) or a mirror
(reflecting top). Across MeshBlock boundaries the inflow intensity is that of the previous sweep
(block Jacobi).

On 3-D spherical-polar and cubed-sphere meshes with one MeshBlock along $x_1$, the rays start
at an adaptive **ray cut** (`clean_ray_cut`, on by default where it applies). Below the cut the
intensities are set to the diffusion limit, whose error is $O(R^2)$ with
$R=|\nabla S|/(\chi S)$ (`CleanRayCut`). The keys `clean_ray_cut_r`, `_margin`, `_wait` and
`_taumin` control it.

### 3.2 Direction sets

| keys | set | notes |
|---|---|---|
| `vet_gd_gl_nmu = N`, `vet_gd_gl_nphi = M`, `vet_gd_gl_stagger` | product set: $N$ Gauss–Legendre nodes in $n_z$ times $M$ uniform azimuths; with `stagger`, alternate rings are shifted by half an azimuth step | weights positive; $\sum w=1$, $\sum w\hat n=0$, $\sum w\,\hat n\hat n=\mathsf I/3$ exactly (`VgdProduct`) |
| `vet_gd_ls = 6, 8, 10, 12` | level-symmetric $LQ_N$ (Carlson), the Athena++/Jiang (2021) family; `vet_gd_ls = 10` has 120 directions | octahedral symmetry |
| `vet_gd_geo = n` | geodesic set: vertices of a subdivided icosahedron ($10n^2+2$ directions), Voronoi solid-angle weights | default of a fresh 3-D run: `vet_gd_geo = 3` (92 directions) |
| `vet_gd_nside = n` | HEALPix ($12n^2$) | the default only on restarts that name no set |

If several set keys are named, the product set takes precedence, then `vet_gd_ls`, then
`vet_gd_geo` (`VgdBuildSet`). A fresh 2-D axisymmetric run defaults to GL $9\times n_\phi$
(staggered), with $n_\phi$ chosen so that $2\pi/n_\phi$ divides the $\phi$ extent. The
reference production decks use GL $7\times18$, staggered (126 directions).

**Rotation (optional).** With `vet_gd_rotate_dt` (> 0 or `auto`) the set is rotated on the
physical clock. `vet_rotate_avg = K` then averages the closure over the last $K$ rotated sets
by carrying their moment offsets (`rad_m1_rotavg.cpp`). `vet_gd_rotate_xfade` cross-fades
between consecutive sets, and `vet_rotate_avg_ema` replaces the ring by an exponential moving
average. A fresh run that rotates defaults to $K=4$. The reference production decks switch
rotation off (`vet_gd_rotate_dt = 0`, `vet_rotate_avg = 0`, `clean_rot_ring = off`).

### 3.3 What the rays hand to the moment solve

Per cell, in the cell's local basis (`CleanClosureFromGd`, `CleanClosureFromRays`):

- $J=\langle I\rangle$, $\mathsf K=\langle I\,\hat n\hat n\rangle$, and
  **$\mathsf D=\mathsf K/J$**, the Eddington tensor (diagonal and off-diagonal parts);
- the half-range flux ratios $h^\pm_b=H^\pm_b/J$, where
  $H^\pm_b=\langle I\,\mu_b\rangle_{\mu_b\gtrless0}$. They satisfy
  $0\le h^+\le1$ and $-1\le h^-\le0$ and equal $\pm\tfrac14$ for an isotropic field;
- $C_1=\langle|I-S|\rangle/J$, a measure of non-locality;
- (staggered consistency term) the rays' own flux through each face, $F_{\rm ray}$, from the
  upwind face intensities (`CleanGdFaceSums`, `CleanFaceHalf`).

`clean_closure = eddington` replaces all of these by the isotropic values
($\mathsf D=\mathsf I/3$, $h^\pm=\pm\tfrac14$, $C_1=0$) without a formal solution. This is useful
for tests.

### 3.4 When the closure is rebuilt

`clean_closure` (default `extrap1`) controls the rebuild (`CleanClosureExtrap`,
`CleanClosureAfterSolve`):

- **extrap1.** In a two-stage step, pass 1 of the stage-1 solve sweeps the rays of its own
  Newton iterate (`CleanClosureIterate`). Stage 2 uses the straight line through the last
  step's and this step's stage-1 moments, evaluated at $t^{n+1}$ and weighted by $1/(1+z)$,
  where $z=c\Delta t/\Delta x$. A cell keeps its latest moments wherever the extrapolated ones
  would not be those of a non-negative intensity. A backward-Euler step instead sweeps the rays
  of its converged solution. So a step normally costs **one ray sweep**.
- `iterate` re-sweeps the rays at every pass, up to `clean_maxpic` (12) passes. `extrap` and
  `extrap2` are older line-extrapolation variants. `predict` builds the closure once per step
  from a predicted $t^{n+1}$ state (`CleanPredictState`).

The closure enters the moment solve through $\mathsf D$ in the stress terms, through
$h^\pm$ and $w(C_1)$ in the kinetic share of the face flux, and through $F_{\rm ray}$ in the
consistency term (§4).

---

## 4. Spatial discretisation (staggered layout)

`clean_layout = staggered` is the only layout; the cell-centred layout was removed. Each cell
holds $E_i$. Each face holds the scalar face momentum $F_f$, the flux along the face's unit
normal. A cell's flux vector is the mean of its two faces in each direction on Cartesian and
spherical-polar meshes, and an exact linear map of its face values on the cubed sphere
(`CleanStagCellVector`, `CleanStagCsVector`). The cell vector is used for output, for the closure
input and for the β terms.

### 4.1 The cell (energy) equation of a stage

For each cell,

$$
E_i\,(1+h c\,k_{E,i}) + \frac{h}{V_i}\sum_{f\in\partial i} A_f\,\Phi_f(E)
\;=\; E^{\rm old}_i + h c\,k_{P,i}\,a T_i^4\,(1-\tfrac12\beta_i^2) + h\,s^\beta_i - W_i ,
$$

where $E^{\rm old}$ is the stage's old vector (§5), $\Phi_f$ is the face **energy flux**, $W_i$
is the work of the radiation force on the gas, and $s^\beta$ is the β source of §2.1. The
exchange is linearised about the Newton iterate: the emission follows $E$ through the gas
temperature, so the $E$ row's net absorption is $(1-r)\,hc\,k_E E$, with $r=s/(\rho c_v+s)$ and
$s=4hc\,k_P aT^3$, and the right-hand side carries the matching $-r\,hc\,k_E E_k$ (§6.3). Spherical-polar and cubed-sphere meshes use
the exact face areas $A_f$ and volumes $V_i$ of the AthenaK coordinates.

### 4.2 The face energy flux

The face momentum equation of a stage,

$$
\frac{F_f-G_f}{h} + c^2\,\partial_n(\mathsf D E)_f = -c\,\chi_f\,(F_f-\mathcal A_f)\;(+\ \beta\ {\rm terms}),
$$

where $G_f$ is the explicit part (old face flux plus stage history) and
$\mathcal A_f=(\mathbf a E)_f\cdot\hat n_f$, can be solved for $F_f$ face by face. This
elimination of $F$ is what turns the stage into a **scalar system for $E$**. The result is
blended with a kinetic (upwind discrete-ordinates) flux taken from the rays:

$$
\Phi_f \;=\; (1-w_f)\Big[\theta_f\big(G_f - h c^2\,\partial_n(\mathsf D E)_f\big) + (1-\theta_f)\,\mathcal A_f\Big]
\;+\; w_f\,c\,\big(h^+_f E_L + h^-_f E_R\big) \;+\;\Phi^{\rm lag}_f ,
\qquad \theta_f=\frac{1}{1+h c\,\chi_f}.
$$

(`CleanFaceWeights`, `CleanFaceSource`, `CleanFaceFlux` in `rad_m1_clean.hpp`.) The three shares:

- **Diffusion / wave share**, weight $1-w_f$. In the thick limit it is the compact diffusion
  flux $-(c/\chi_f)\partial_n(\mathsf D E)$. In the thin limit it is the staggered wave operator.
  It needs no Rhie–Chow term or added dissipation. On spherical-polar meshes the normal stress
  carries the curvature terms ($2P_{rr}/r$ and related; `CleanStagSph`, `CleanFaceMetric`).
- **Kinetic share**, weight $w_f$: the first-order upwind discrete-ordinates flux in moment
  form. $h^+_f\ge0\ge h^-_f$ are the rays' half-range ratios reconstructed to the face from the
  upwind side (minmod-limited). With `clean_face_jse` (default on) they take the shape of the
  cells' half-range fluxes at the face (`CleanFaceCoefficients`). This share is an M-matrix in
  any dimension.
- **Lagged terms** $\Phi^{\rm lag}_f$: the consistency term, the off-diagonal stress and the
  deferred part of the advective flux (below). All of them are built at the previous iterate,
  so the operator stays a scalar M-matrix.

**Thin-share weight.** With `clean_w_fse` (default on),

$$
w=\mathrm{smoothstep}\!\Big(\frac{g-c_0}{c_1-c_0}\Big),\qquad g=\max\!\big(|H|/J,\ \epsilon\,C_1/1.5\big),
$$

with $c_0,c_1$ = `clean_w_c0`, `clean_w_c1` (0.3, 0.6) and $\epsilon=\min(k_P/\chi,1)$. With
`clean_w_fse_tau` (default on), scattering isotropises only over the cell's optical depth
$\tau_d=\chi V/A_f$, so $\epsilon\to1-(1-\epsilon)(1-e^{-\tau_d})$. The kinetic share therefore
switches on where the field is anisotropic, not merely where it is non-local.

**Face opacity.** $\chi_f$ is the mean of the two cells. With `clean_face_opac = harmonic`
(default) an $x_1$ face takes the distance-weighted harmonic mean, so the transparent side
limits the flux at an unresolved opacity jump (`CleanFaceChi`).

### 4.3 Consistency term and `clean_face_hctau`

`clean_consistency` (default on wherever the rays exist) adds a lagged correction that makes
the scheme's thin-limit flux match the rays' own transport solution. Two forms exist:

- `clean_stag_thinf = false` (the reference production form): a lagged additive flux on the
  kinetic share,
  $\Phi^{\rm lag}_f\mathrel{+}=w_f\,q_f\,c\,(F_{y,f}-h^+_fJ_L-h^-_fJ_R)$ (`CleanThinFaceDG`).
  Here $F_{y,f}$ is the rays' face flux clamped to $[-J_R,J_L]$, and $q_f$ is a gate.
- `clean_stag_thinf = true` (the code's fresh default): the same remainder is folded into the
  face half ranges, $g^\pm$, kept within $[0,1]$ and $[-1,0]$. This form is realisable by
  construction.

With `clean_face_hctau = true` (on in production; default off), the face value of the rays'
flux is interpolated linearly in optical depth between the two cell centres:

$$
H_f=\frac{\tau_R H_L+\tau_L H_R}{\tau_L+\tau_R},
$$

where $\tau$ are the half-cell depths, and $H$ is scaled by $r^2$ on radial faces of curved
meshes (exact for free streaming). Next to an opaque cell the face then sees the flux the
surface emits. The rays' own face flux is kept where the two cells' $H$ jump, for instance at a
beam edge or a shadow (`CleanHcTauW`, `CleanStagBracket`).

### 4.4 Advective flux

The share $(1-w)(1-\theta)\mathcal A_f$ is central by default: $a_f=(a_L+a_R)/2$. On a face
where the central form would break the M-matrix (cell Péclet number above 1), the smallest
donor-cell share $\varphi$ that restores it is taken implicitly. The remainder
$\varphi\,[\,\text{high-order}-\text{donor}\,]$ is carried as a deferred correction at the
iterate $E^k$ (`CleanStagAdvectCorr`). The high-order remainder has three forms:

| key | face value carried at the fixed point |
|---|---|
| (neither key) | $a_f E_f$ with $E_f$ the van Leer (PLM) value from the upwind side |
| `clean_stag_upwind_central = true` | $a_f E_{\rm up}+a_{\rm down}\,m/2$, with $m={\rm minmod}(2(E_{\rm up}-E_{\rm upup}),\,E_{\rm down}-E_{\rm up})$: the central mean of the cells' $aE$ within Sweby's TVD bound |
| `clean_stag_central = true` (implies the above) | the plain central mean $\tfrac12(a_LE_L+a_RE_R)$. The radiation work's drain enters the $E$ rows' diagonal, and the gas kick enters with the same Patankar weight (`CleanStagMprkWeight`), so positivity is kept. It needs `clean_xt_exact = cons`. |

Both keys are on in the reference production decks; their code default is off. The reason, from
the code comments: the van Leer remainder is upwind-biased at $O(\Delta x^2)$, and in convective
zones its up- and downflow parts do not cancel.

### 4.5 Off-diagonal stress (`clean_stag_xdonor`, `clean_stag_xsrc`)

The transverse stress flux $C=-(1-w)\theta h c^2\,\partial_t(D_{nt}E)$ would give a 9-point
stencil that is not an M-matrix. With `clean_stag_xdonor = lag` (default), $C$ is evaluated at
the iterate and entered as a donor-cell coefficient $q\,E_{\rm donor}$ with
$q=C(E^k)/E^k_{\rm donor}$ and $|q|\le c$. The part whose explicit outflow over the stage stays
below $\tfrac16E^{\rm old}_{\rm donor}$ stays explicit, and the rest becomes donor-implicit
(`CleanStagXDonor`, `CleanStagXSplit`). At the fixed point the face carries the central $C$
wherever $|C|\le cE_{\rm donor}$. `clean_stag_xsrc` (default on) applies the same treatment to
the face source terms. From pass `clean_stag_xdonor_freeze` (default 4) on, the donor side and
the cap are held, which removes period-2 flip-flops of the passes. The gas momentum always uses
the central stress.

### 4.6 The face momentum $F_f$ (the F solve)

Once $E$ of a pass is known, $F_f$ is solved from the face momentum equation with E fixed,
one unknown per face. A dissipation term is added on the stress:

$$
F_f\Big(\frac1{\theta_f}+\frac{hc(\alpha_L+\alpha_R)}{2\Delta x_f}\Big)-\frac{hc\,(\alpha_LF_{f-1}+\alpha_RF_{f+1})}{2\Delta x_f}=\frac{F^H_f}{\theta_f},
$$

where $F^H$ is the eliminated flux of §4.2 and the weight
$\alpha=\max(\mathcal A(\tau)^2,\ \ldots)$ with $\mathcal A(\tau)=(1-e^{-\tau})/\tau$
(`CleanDissipationWeight`). The dissipation is therefore $O(1/\tau^2)$ in thick cells. On the
streaming share, the kinetic upwind flux of the face momentum replaces the cell momentum flux,
so a beam keeps $r^2F$. The lines along each direction are tridiagonal M-matrices. They are
solved by Thomas or parallel cyclic reduction (`clean_stag_fs_line = pcr | thomas`).
`clean_stag_fschur` (default where it applies) joins the lines across MeshBlocks and periodic
wraps exactly. (`CleanStagFSolve`, `CleanStagFRow`, `CleanStagFSchur`.)

Three options bound or balance $F$:

- **`clean_stag_fwb`** (default on). The F rows are written in deviation form about the rays'
  own state $F_{\rm ray}$ (the rays' flux factor times the donor $E$, plus the advective part).
  $F_{\rm ray}$ is then the *exact* steady state of the streaming rows, whose stress, curvature
  and transport would otherwise differ from the rays' through discretisation. The residual is
  weighted by $1/(1+\tau_f)$ and is radiation momentum: it is recorded as such and not exchanged
  with the gas (§7.2). `clean_stag_fwb_stage`, `_t2`, `_line` and `_ring` choose when the
  residual is formed and project it (§7.3).
- **`clean_stag_fcone`** (default on with fwb). Each cell's flux vector is kept in its
  realisable cone, $|\mathbf F_i|\le cE_i$, by scaling faces towards 0 with
  $l_f=\min(l_L,l_R)$ and $l_i=\min(1,cE_i/\sum_f|F_f||u_f|)$ (`CleanStagFCone`). The cut
  radiation momentum is recorded and printed.
- `clean_stag_flim` (default off): an older bound on $F^H$ by the donor's $cE$. It also cuts
  free streaming.

`clean_stag_vimp` (default on; on the cubed sphere also `clean_stag_vimp_cs`) makes the gas
velocity in the face momentum implicit through a drag share (`CleanStagDrag`).

### 4.7 Boundaries

Allowed `implicit_bc_*`: `marshak`, `outflow`, `reflect`, `periodic`, and on $x_1$ an imposed
`flux`. With the default `clean_bc = halfrange`, an open $x_1$ end carries the end cell's
half-range flux $c\,h^\pm E$. `clean_bc_area` (default on) writes it in area form, so an isotropic
field is balanced by the bath. A wall carries no flux, and an imposed flux is the face's whole
energy flux and momentum (`CleanFaceWeights`, `CleanFaceSource`). The vet_gd sweep needs
`implicit_bc_x1max = marshak` or `reflect`; transverse Marshak faces are not supported.

---

## 5. Time integration

### 5.1 IMEX structure

`<rad_m1>/time_scheme = hesdirk2` (the default where accepted) needs `<time>/integrator = rk2`.
The explicit part *is* AthenaK's Heun hydro, unchanged. Each implicit stage is a
backward-Euler-type solve of length $h=\gamma\Delta t$, $\gamma=1-1/\sqrt2$, on a supplied old
vector (`rad_m1_time2.cpp`; tableau in `src/utils/imex_sdirk2.hpp`):

| `time2_tableau` | implicit tableau | stage-1 old vector | stage-2 old vector |
|---|---|---|---|
| `hesdirk2` (default) | ESDIRK, FSAL: $A=[[0,0,0],[1-\gamma,\gamma,0],[\tfrac12,b_2,\gamma]]$, $b_2=\gamma/(2(1-\gamma))$ | stage start $+(1-\gamma)\Delta t K_1$ | Heun average $+\Delta t(\tfrac\gamma2K_1+(b_2-\tfrac\gamma2)K_2)$ |
| **`sdirk2`** (reference production) | Alexander's two-stage L-stable SDIRK, $b=(1-\gamma,\gamma)$; no explicit $K_1$ | gas at $t^n+\gamma\Delta t$ (below) | Heun average $+\Delta t\,[(1-\gamma)-\tfrac\gamma2]\,K_a$ |
| `trbdf2` | ESDIRK form of TR-BDF2 | (M1 tests only) | |

A stage slope is $K=(Y-\text{old vector})/h$.

**sdirk2 with the Heun hydro** (`Time2GasTime`). Let $U^*=U^n+\Delta t\,L(U^n)$ be Heun stage 1.
Stage 1 solves on $S_1=U^n+\gamma(U^*-U^n)$, the gas at $t^n+\gamma\Delta t$. This choice gives
each stage one time for both parts, so the IMEX coupling conditions
$b_E\cdot c_I=b_I\cdot c_E=\tfrac12$ hold. After the solve the gas becomes $U^*+d$ with
$d=\gamma\Delta t\,K_a$. Heun stage 2 evaluates its hydro fluxes on the "flux state"
$U^*+f\,d$, where

$$
f=\frac{1+\gamma z}{\gamma(1+z)},\qquad z=\Delta t\Big(c\,k_E+\frac{4c\,k_P aT^3}{\rho c_v}\Big).
$$

So $f=1/\gamma$ where the exchange is not stiff and $f\to1$ where it is. With `time2_gtime_mom`
(default on) the momentum always takes $f=1/\gamma$, because the radiation force is not a
relaxation of $v$. The flux state's internal energy follows its kinetic-energy change. The
stage-2 implicit solve then starts from the Heun average plus $[(1-\gamma)-\gamma/2]\Delta t K_a$.
Together these give the total implicit weights $1-\gamma$ on $K_a$ and $\gamma$ on $K_b$.

### 5.2 Order of operations in one step (sdirk2)

1. Save the stage-start state; `dt_retry` holds a full cycle checkpoint.
2. Heun hydro stage 1 (including gravity and sources) gives $U^*$.
3. Set the gas to $S_1$ (`Time2GasTime(1)`) and form the stage-1 old vector (`Time2FormStage`).
4. Radiation stage 1 (`m1_stagen` task list): closure limits, opacities, then `CleanSolve`
   (§6), then primitive recovery and ghost exchange of the moments.
5. Build the hydro flux state $U^*+f\,d$ (`Time2GasTime(3)`).
6. Heun hydro stage 2 (Heun average).
7. Radiation stage 2: `CleanSolve` on the stage-2 old vector. The result is $U^{n+1}$.
8. Closure limits, the new radiation time step, and with `time2_cfl_post` the hydro CFL of the
   coupled state.

These lists are separate task lists (`rad_m1_tasks.cpp`) that the driver runs after each hydro
stage (`src/driver/driver.cpp`). A run without `<rad_m1>` is unchanged.

### 5.3 Fallbacks

- **Backward Euler.** The step is taken as backward Euler (one solve with $h=\Delta t$, storing
  $K_1$) when no valid slope exists: the first step (`time2_start = be`), a restart from a file
  without the slope, or the step after a failed stage. With `time2_auto_be = N` (default 16)
  the next 1, 2, 4, …, N steps after a stage fallback go straight to the backward-Euler path.
- **Stage failure.** A stage fails when its solve did not converge, left $E<0$ that the
  conservative floor did not absorb, or left gas $e<0$. The cycle is then put back and redone by
  backward Euler (`Time2Restore`). With `time2_reject = retry` it is retried at a smaller $\Delta t$
  instead, and after `<time>/dt_retry_t2_be` such cuts by backward Euler.
- **`<time>/dt_retry`** (default on for `<hydro>`/`<mhd>` + `<rad_m1>` on a uniform mesh). A
  failed backward-Euler solve, or a write-back that leaves $e<0$, restores the cycle checkpoint
  and retries with $\Delta t\times$`dt_retry_fac` (0.5; `dt_retry_fac_t2` for stage
  rejections). $\Delta t$ grows again by `dt_retry_grow` after `dt_retry_ngood` good cycles.
  Without `dt_retry`, a non-converged backward-Euler solve is fatal if its last change exceeds
  `implicit_resid_fatal`.

### 5.4 Stage options

| key | default | effect |
|---|---|---|
| `clean_predict` / `implicit_predictor = step` | on | first Newton iterate $E^{\rm old}+hK_1$ where that changes $T_{\rm rad}$ by less than `implicit_predictor_tcap` (`CleanFirstIterate`); affects the starting point only |
| `clean_pred_stage` | off | stage 2 starts from $E^{\rm old}+hK$ of *this* step's stage 1 instead of the last step's stage 2 |
| `clean_stag_eicarry` | off (on in production) | a negative gas internal energy at the stage start enters the gas row as it is, paid by the implicit exchange, instead of being reset to 0 (which would create energy) (`CleanOldVector`) |
| `time2_k1_stiff = s` | 0 (off) | where $k\Delta t>s$, moves the exchange part of $K_1$ from the gas to $E$ (energy-neutral), so the stiff exchange is taken by the implicit stages only. Relevant to tableaus that use $K_1$ (`hesdirk2`). |
| `time2_face_fct` | on (fresh clean run) | Zalesak limiter on the face part of the explicit stage increment: a cell whose $E^{\rm old}$ would not stay positive has its outflowing faces scaled by its donor ratio $R$ on both sides (energy exact). It keeps `time2_face_fct_emin` (0.1) of $V E^n$ in a limited cell. |
| `time2_enth_be` | on (fresh) | backward-Euler solves use the stage rule for the advective face coefficient |

### 5.5 Steady start

`clean_init = steady` relaxes the initial radiation to the initial matter with the gas frozen.
It takes backward-Euler radiation steps of `clean_init_dtfac`×Δt (default $10^6$) until
$\max(|\Delta E|,|\Delta F|/c)<$ `clean_init_tol` (1e-8) of max $E$, at most `clean_init_n`
passes. A step that fails is undone and shrunk. Because $\Delta t\to\infty$ drops the time term,
the result is the scheme's own discrete steady state, including rays, closure and consistency
terms (`CleanSteadyStart`). The closure history is then dropped.

---

## 6. The implicit stage solve (`CleanSolve`)

### 6.1 Structure

```
CleanSolve(stage):
  closure for the stage (extrapolated, §3.4); old vector (stage increments added);
  face coefficients (w, theta, G, h+-); consistency term
  for pass = 0 .. clean_maxnl-1:
      [stage 1, pass 1: closure of this iterate (one ray sweep)]
      [pass < clean_stag_opac_passes: opacities at the iterate's T]
      xdonor coefficients q of the iterate (held from clean_stag_xdonor_freeze on)
      right-hand side, exchange linearised about (E_k, T_k)
      preconditioner set-up
      X = BiCGStab(A, B, x0 = E_k)                       -- scalar E system
      Newton update: E_k = X; T_k = root of the gas row   -- per cell
      velocity update: F solve (lines), gas momentum, advective coefficients
      convergence test (budget norm)
  realisability count (E < 0; face flux beyond c E_donor)
  write-back: face energy fluxes (Zalesak-limited), E, F, gas momentum and energy, slopes K
```

### 6.2 Linear solver

The operator $A$ is the left side of §4.1 with the eliminated fluxes. It is a non-symmetric,
column-dominant M-matrix. It is solved by **right-preconditioned BiCGStab** (`CleanBiCGStab`,
`rad_m1_clean_solver.cpp`), with breakdown restarts and a check on the true residual. The
preconditioner (`clean_stag_pc`, default `mg`):

- `line`: an exact tridiagonal solve along every $x_1$ line, plus the diagonal of the
  $x_2$/$x_3$ couplings;
- `mg` (default): a multigrid V-cycle. Its smoother is line relaxation over the $x_1$ lines
  (PCR in team scratch on GPUs, Thomas on the host), with the lines coloured either
  mirror-symmetrically about the MeshBlock's mid-planes (`clean_stag_mg_color = mirror`,
  default) or by $(j+k)$ parity. The smoothing sweep is set by `clean_pc_sweep` (`sym` on a
  fresh run, `rb` in the reference decks). Coarse levels use Galerkin $RAP$ on $2\times2$
  $(x_2,x_3)$ aggregates within a MeshBlock, with volume weights on curved meshes
  (`CleanStagVCycle`, `CleanStagCoarseSetup`). The preconditioner is MeshBlock-local, unless
  `clean_stag_mg_tile` is set (§7.8).

The operator is applied from stored stencil rows (`clean_stag_stencil`, default on). Halo
exchange overlaps the inner cells.

### 6.3 Gas coupling: Newton in E, exact root in T

In pass $k$ the emission is linearised about $T_k$: its response to $E$ enters the $E$ row as
$(1-r)\,hc\,k_E E$, with $r=s/(\rho c_v+s)$ and $s=4hc\,k_P aT_k^3$. After the linear solve each
cell's gas temperature is the unique root of its local energy row,

$$
e(\rho,T) + hc\,k_P\,(1-\tfrac12\beta^2)\,a\,T^4 \;=\; e^{\rm old} + hc\,k_E\,E_{k+1} + (\text{work, β terms}),
$$

found by safeguarded Newton with bisection fallback for a general EOS (`CleanExchangeTemperature`,
tolerance `clean_troot_rtol` = 1e-13). With `clean_xt_exact = cons` (default) the first
iterate's $T$ comes from the pass's own gas row. The work and kinetic-energy mismatch of the
linearisation is booked to $E$, so the gas energy takes the heating only and **total energy closes
exactly** (`CleanRightHandSide`, `CleanNewtonUpdate`). If a solve's change stalls (minimum over
the last `clean_gas_newton_switch` = 4 passes not below `clean_stall_fac` = 0.5 of the earlier
minimum, rising at least twice, from pass `clean_gas_newton_switch_min` = 6), the remaining
passes drop the Newton slope and iterate the bracketed root as a contraction.

### 6.4 Convergence

With `clean_conv = budget` (default), every row is measured in the cell's own energy budget for
the stage, $S_i$, the sum of its gross terms: old value, face fluxes, exchange, work and β terms.
A change is therefore a fraction of the energy the cell exchanges. The nonlinear tolerance is

$$
{\rm tol} = \min\!\big(\max(\texttt{clean\_conv\_kappa}\times\varepsilon_t,\ \texttt{clean\_tol\_nl}),\ \texttt{clean\_conv\_max}\big),
$$

where $\varepsilon_t=\tfrac h2|f(Y)-f(U^n)|$ is the step's time-error estimate in the same weights.
The floor is at least $10\times$`clean_tol`. BiCGStab runs in the same norm with
**Eisenstat–Walker forcing** (choice 2, safeguarded):
$\eta=\min(0.1,\ 0.9\,(r_k/r_{k-1})^2)$, with $\eta=0.1$ on the first pass. ($r_k$ is the
pass's initial residual. The code's safeguard $\max(\eta, 0.9\eta_{k-1}^2)$ binds only when
the cap 0.1 does.) The linear tolerance is
$\max(\eta r_0,\ 0.1\,{\rm tol})$. The forcing constants are fixed in the code;
`implicit_lin_ew_max` is not read by this solver. A solve that is not in a budget pass (a single
solve, or the first cold step) uses `clean_tol` relative to $|B|$ together with the per-row test
`clean_lin_cnorm`. `clean_conv = strict` measures $E$ and $T$ relative to themselves against
`clean_tol_nl`.

| key | default | meaning |
|---|---|---|
| `clean_tol` | 1e-10 | linear relative tolerance (strict passes) |
| `clean_tol_nl` | 1e-11 | nonlinear floor (raised to ≥ 10 `clean_tol`) |
| `clean_conv_max` | 1e-5 | ceiling of the budget tolerance |
| `clean_conv_kappa` | 1.0 | factor on the time-error estimate |
| `clean_maxnl` | 200 | max passes per solve |
| `clean_maxit` | 800 (production 200) | max BiCGStab iterations per pass |
| `clean_lin_pass` | on (fresh) | a linear solve that hits `clean_maxit` hands its iterate to the next pass |
| `clean_lin_cnorm`, `clean_init_cnorm` | 1e-10 | per-row residual test $\max_i|r_i|/(a_{ii}E_i)$ |

---

## 7. Conservation and robustness properties

The properties below are design guarantees of the code paths named. They are exact to
round-off unless stated otherwise.

### 7.1 Energy (gas + radiation)

- $E$ is updated in flux form with one energy flux $\Phi_f$ per face, shared by both cells.
  MeshBlock-shared faces take their owner's value (`CleanStagFaceTie`).
- The gas takes exactly the exchange and work that $E$ gives up (`CleanWriteBack`;
  `clean_xt_exact = cons`, §6.3).
- Both Zalesak limiters (`time2_face_fct` and `clean_wb_fct`) scale a face's flux on both sides
  at once.
- The energy that crosses the physical $x_1$ boundaries is accumulated for the history output
  (`CleanBndAcc`).

Conservation is broken only by explicit floors: `e_floor` raises, the hydro's density and
pressure floors, and the `floorcap_*` caps. These are counted and reported, and in the reference
configuration `clean_wb_fct` makes $E$-floor raises rare.

### 7.2 Momentum

The gas receives, face by face, minus the radiation momentum change across the face (§2.2). The
stress term telescopes, so total momentum changes only through boundaries. There are two
documented exceptions. Both act on the streaming share and both are recorded as radiation
momentum, not exchanged with the gas: the well-balancing residual of `clean_stag_fwb`, and the
realisability cut of `clean_stag_fcone` (or `flim`). Where they act, gas + radiation momentum is
not conserved exactly; the size of the imbalance is reported (`CLEAN FCONE` lines and the
ledgers).

### 7.3 Angular momentum $L_z$

- **3-D spherical polar.** The lateral $\phi$-face stress is written in torque form
  (`CleanFaceStress`). `clean_stag_fwb_ring` projects the fwb residual of each $\phi$ ring onto
  zero net torque. It is the code default on a uniform 3-D mesh with periodic $x_3$, but the
  reference production decks set it to 0.
- **Cubed sphere.** $L_z$ is kept by telescoping (`clean_cs_lz_tel`, `CleanCsTorque`), and the
  gas $\phi$ momentum takes the faces' absorption.
- **Ledger.** Whatever $L_z$ the radiation keeps (fwb or fcone residuals, radial-boundary
  torque) goes into a ledger stored in the restart parameters (`CleanLzBook`). Gas +
  radiation + ledger (history column `Lz_tot`) is conserved to round-off on closed test
  problems (2e-16 over 40 cycles, `tests/m1_clean` `conservation_*`). Gas + radiation alone
  drifts by what the ledger holds: 7e-7 on the spherical-polar shell (the `clean_stag_fcone`
  cut; 2e-15 without it) and 6e-5 on the cubed-sphere shell. So $L_z$ is exact by
  construction only where no cut is booked.

### 7.4 Mass

The radiation solve never changes $\rho$. Mass conservation is that of the hydro scheme and its
floors.

### 7.5 Positivity and realisability

- The $E$ operator is an M-matrix for any $\mathsf D$ and $h$. This rests on four constructions:
  the kinetic share with $h^+\ge0\ge h^-$, the minimal donor share on the advective flux, the
  donor-cell form of the lagged off-diagonal stress (`clean_stag_xdonor`), and the implicit
  upwind memory (`CleanStagMemUp`). A non-negative right-hand side therefore gives $E\ge0$.
  The explicit parts of the right-hand side are bounded, so the right-hand side stays
  non-negative.
- `time2_face_fct` limits the explicit stage increment, and `clean_wb_fct` (default on, fresh)
  limits the written energy fluxes, so that $E>e_{\rm floor}$ by construction (`CleanWbFct`).
- `clean_stag_fcone` keeps $|\mathbf F_i|\le cE_i$ by construction.
- Two realisability counts are reported each solve (`CleanStagRealizable`): R1, cells with
  $E<0$, and R2, interior faces whose flux exceeds $cE$ of their donor.
- A write-back that leaves gas $e<0$ is treated as a failed solve (backward-Euler redo or
  `dt_retry`). It is not clipped silently.

### 7.6 Well-balancing

Two mechanisms work together:

1. **Hydrostatic + radiative equilibrium at rest.** With
   `force_reference = wb_arad`, the problem generator supplies a reference radiative
   acceleration $\mathbf a_{\rm ref}$. The hydro's well-balanced scheme already includes it in
   the gas momentum, so the coupling hands the gas only the *residual* force
   $\Delta(\rho\mathbf v)-h\rho\mathbf a_{\rm ref}$. The reference work $h\rho\,\mathbf a_{\rm ref}\cdot\mathbf v$
   is paid by $E$ (`force_reference_work = split`, which the clean scheme requires). Other
   choices: `clean`, which takes the clean solve's own force from its first coupled solve, and
   `wb_arad_md`, which captures the steady start's face forces in multi-D (needs
   `clean_init = steady`).
2. **Radiation steady state.** `clean_stag_fwb` makes the rays' flux the exact discrete steady
   state of the streaming face rows (§4.6). The consistency term (§4.3) and the steady start
   (§5.5) make the scheme's own steady solution the initial state.

### 7.7 Restart and `dt_retry`: bitwise

- **What is stored.** The restart file carries the full scheme state: face momenta and
  energy-flux memory, the stage slope $K_1$, closure moments for the extrapolation, the ray
  cut, ring offsets and the $L_z$ ledger (`Time2RstPack`, format version 7). A run of $N$ steps
  is designed to be bitwise identical to $N/2$ steps, a restart and $N/2$ more steps. The
  development gates test this property.
- **Keys and defaults.** Every key a run uses is recorded in the restart's parameter dump.
  **A restart takes the stored value, not the code's current default.** A restart from a file
  that lacks a key takes that key's *old* default, and the log echoes it. Name every key that
  matters in the restart input.
- **`dt_retry`.** The cycle checkpoint (`src/driver/dt_retry.cpp`, `dt_retry_views = minimal`)
  restores every evolving array, so a retried cycle is the cycle at the cut $\Delta t$.
- **Determinism.** The solver's decisions (pass counts, stall switches, limiter rounds) are
  functions of the solve's own history and of globally reduced quantities. Non-finite
  reductions are mapped to $+\infty$ on every rank (`CleanMaxRanks`).

### 7.8 MeshBlock-layout and rank-count invariance

`<mesh>/layout_invariant` (default true on fresh runs) routes every global or shell sum that
feeds the state through an exact binned accumulator (`src/utils/repro_sum.hpp`,
`rad_m1_layreduce.hpp`). Such a sum depends only on the multiset of its terms. It also computes
spherical-polar $\theta,\phi$ from the global cell index.

| `<mesh>/layout_invariant_mode` | sums | what is guaranteed | cost |
|---|---|---|---|
| `block` (default) | each MeshBlock's partial in a fixed order, then the partials summed exactly | independent of the rank count and the MeshBlock→rank map; still depends on the MeshBlock size, through the partials and the MeshBlock-local preconditioner, at the solver-tolerance level | ≈ 0 |
| `exact` | every cell term deposited exactly | sums independent of the layout; the preconditioner is still MeshBlock-local | ≈ +8 % per cycle (measured on one production case) |
| `tile` + `<mesh>/layout_tile = N` | block mode with fixed $N\times N$ $(x_2,x_3)$ tiles as the summation units; history and shell sums exact; sets `clean_stag_mg_tile = N`, so the multigrid acts on tiles | sums **and** preconditioner identical for every layout whose MeshBlocks are whole tiles; needs one MeshBlock along $x_1$ and $N$ dividing the MeshBlock's $n_{x2}$ ($n_{x3}$) | ≈ block |

The exactly joined cross-MeshBlock F line solves (`clean_stag_fschur`) still differ at
round-off between layouts in this release (see §12): on the test shell in tile mode, one
MeshBlock vs 1 x 2 (2 x 2) MeshBlocks differ by 1.5e-10 (4.5e-10) of the largest value after 8
cycles (`tests/m1_clean` `layout_tile`).

---

## 8. Supported geometries

| mesh | closure | notes |
|---|---|---|
| Cartesian (1-D, 2-D, 3-D) | `vet_sc` short characteristics, or `clean_closure = eddington` | uniform cells; `clean_consistency` on by default |
| spherical polar | `closure = vet_col` with `vet_gd = true` (the fresh default) and `clean_rays = sc` | 3-D wedge, 2-D axisymmetric (the $\phi$ extent must map the direction set onto itself), or full sphere with `use_polar_boundary`; radial stretching allowed; several MeshBlocks along $x_1$ only on a wedge; ray cut only with one MeshBlock along $x_1$ |
| cubed sphere | as spherical polar; `rad_m1_vetgd_cs.cpp` | $F$ along each face family's normal; $L_z$ by telescoping; `clean_cs_lz_tel`, `clean_cs_lz_bnd_free` |

Not supported (FATAL):

- mesh refinement (SMR/AMR) with `scheme = clean`;
- $\hat c<c$;
- `time2_tableau = sdirk2` with MHD;
- transverse Marshak boundaries in the vet_gd sweep.

The vet_col keyword is only the interface that the vet_gd closure reads. The per-column 1-D
solver it once named is not in this build.

---

## 9. Reference production configuration

The `<rad_m1>` keys of the reference stellar-envelope production decks (3-D spherical-polar
wedges, tabulated opacities, a general EOS) that matter for the scheme. Everything not listed
is at its fresh-run code default.

| key | value | code default (fresh) |
|---|---|---|
| `transport`, `scheme`, `clean_layout` | `implicit`, `clean`, `staggered` | same |
| `closure`, `vet_gd`, `vet_gd_replace`, `clean_rays` | `vet_col`, true, true, `sc` | same on sp/cs |
| `vet_gd_gl_nmu`, `_nphi`, `_stagger` | 7, 18, true | `vet_gd_geo = 3` |
| rotation (`vet_gd_rotate_dt`, `vet_rotate_avg`, `clean_rot_ring`) | 0, 0, off | rotating, 4 sets |
| `vet_scatter`, `vet_scatter_form` | true, `ma` | off |
| `time_scheme`, `time2_tableau` | `hesdirk2`, `sdirk2` | `hesdirk2`, `hesdirk2` |
| `time2_gtime_mom`, `time2_face_fct`, `time2_k1_stiff` | true, true, 0 | same |
| `clean_closure` | `extrap1` | same |
| `clean_init` | `steady` (`clean_init_dtfac = 1e6`) | given |
| `clean_stag_central`, `clean_stag_upwind_central` | true, true | false, false |
| `clean_face_hctau`, `clean_face_jse` | true, true | false, true |
| `clean_stag_thinf`, `_thinf_held`, `clean_stag_flim` | false, false, false | true, false, false |
| `clean_stag_fwb`, `clean_stag_fcone` | true, true | same |
| `clean_stag_fwb_stage`, `_fwb_t2`, `_fwb_ring` | false, false, 0 | false, (inert), **true** on 3-D periodic sp |
| `clean_stag_xdonor`, `clean_stag_xsrc`, `_xdonor_freeze` | `lag`, true, 4 | same |
| `clean_stag_eicarry`, `clean_pred_stage`, `clean_stag_kpat` | true, false, false | false, false, false |
| `clean_xt_exact`, `clean_wb_fct` | `cons`, true | same |
| `clean_stag_opac_update`, `_opac_passes`, `clean_face_opac` | true, 3, `harmonic` | same |
| `clean_stag_lat_sums`, `clean_stag_opac_drag`, `clean_stag_tie_gid` | 0, 0, 0 | **true** on a fresh uniform run |
| `clean_stag_mg_color`, `clean_pc_sweep`, `clean_stag_fs_line` | `mirror`, `rb`, `pcr` | `mirror`, `sym`, `pcr` |
| `clean_conv`, `clean_conv_max`, `clean_maxit` | `budget`, 1e-5, 200 | `budget`, 1e-5, 800 |
| `force_reference`, `force_reference_work` | `wb_arad`, `split` | `none` |
| `<time>/integrator`, `cfl_number`, `dt_retry` | `rk2`, 0.3, true | |

---

## 10. Code map

All paths are relative to `src/`.

| concept | file | main function(s) |
|---|---|---|
| key parsing, defaults, FATAL checks | `rad_m1/rad_m1.cpp` | `RadiationM1::RadiationM1` (clean block) |
| clean preflight, allocation, restrictions | `rad_m1/rad_m1_clean.cpp` | `CleanPreflight`, `CleanInit` |
| task lists (stage order) | `rad_m1/rad_m1_tasks.cpp`, `driver/driver.cpp` | `AssembleRadM1Tasks`, stage loop |
| stage solve driver | `rad_m1/rad_m1_clean.cpp` | `CleanSolve` |
| old vector, first iterate | `rad_m1/rad_m1_clean.cpp` | `CleanOldVector`, `CleanFirstIterate` |
| right-hand side, β terms, work | `rad_m1/rad_m1_clean.cpp`, `rad_m1_clean.hpp` | `CleanRightHandSide`, `CleanStagBetaE` |
| face energy flux | `rad_m1/rad_m1_clean.hpp` | `CleanFaceWeights`, `CleanFaceSource`, `CleanFaceFlux`, `CleanFaceChi` |
| face coefficients $w,\theta,h^\pm$ | `rad_m1/rad_m1_clean_solver.cpp` | `CleanFaceCoefficients` |
| consistency term, hctau | `rad_m1/rad_m1_clean.cpp`, `.hpp` | `CleanStagBracket`, `CleanHcTauW`, `CleanThinFaceDG` |
| advective flux (upwind/central) | `rad_m1/rad_m1_clean.hpp` | `CleanStagAdvectCorr` |
| off-diagonal stress | `rad_m1/rad_m1_clean.cpp`, `.hpp` | `CleanStagXDonor`, `CleanStagXSplit`, `CleanFaceCross` |
| F solve, fwb, fcone | `rad_m1/rad_m1_clean.cpp` | `CleanStagFSolve`, `CleanStagFSchur`, `CleanStagFCone`, `CleanStagFwbProject` |
| gas momentum exchange | `rad_m1/rad_m1_clean.hpp`, `.cpp` | `CleanMomentumChange`, `CleanVelocityUpdate`, `CleanStagDrag` |
| gas temperature / Newton | `rad_m1/rad_m1_clean.cpp`, `.hpp` | `CleanNewtonUpdate`, `CleanExchangeTemperature` |
| opacities | `rad_m1/rad_m1_opacity.hpp`, `rad_m1_clean.cpp` | `M1TableOpacities`, `CleanStageOpacity` |
| linear solver, preconditioner | `rad_m1/rad_m1_clean_solver.cpp` | `CleanBiCGStab`, `CleanApply`, `CleanPreconditionerSetup`, `CleanStagVCycle`, `CleanStagLineSolve` |
| write-back, Zalesak on fluxes | `rad_m1/rad_m1_clean.cpp` | `CleanWriteBack`, `CleanWbFct` |
| realisability counts | `rad_m1/rad_m1_clean.cpp` | `CleanStagRealizable`, `CleanCountNegative` |
| L_z ledger | `rad_m1/rad_m1_clean.cpp` | `CleanLzBook`, `CleanCsTorque` |
| closure build / extrapolation | `rad_m1/rad_m1_clean_closure.cpp` | `CleanClosure`, `CleanClosureExtrap`, `CleanClosureIterate`, `CleanClosureAfterSolve`, `CleanRays` |
| ray moments → closure | `rad_m1/rad_m1_clean_closure.cpp` | `CleanClosureFromGd`, `CleanGdFaceSums`, `CleanFaceHalf`, `CleanRayCut` |
| vet_gd sweep (sp), direction sets | `rad_m1/rad_m1_vetgd.cpp`, `rad_m1_vetgd_*.inc` | `VetGdBuild`, `VgdBuildSet`, `VgdProduct`, `VgdLevelSym`, `VgdGeodesic` |
| vet_gd on the cubed sphere | `rad_m1/rad_m1_vetgd_cs.cpp` | (cs band, reads) |
| Cartesian short characteristics | `rad_m1/rad_m1_vet.cpp` | `VetShortChar` |
| ray source function | `rad_m1/rad_m1_vetlat.cpp` | `VetLatSweep` (stage 0) |
| rotation averaging | `rad_m1/rad_m1_rotavg.cpp`, `rad_m1_ringavg.cpp` | `RotAvgParse`, ring moments |
| IMEX stages, tableau | `rad_m1/rad_m1_time2.cpp`, `utils/imex_sdirk2.hpp` | `Time2Init`, `Time2FormStage`, `Time2GasTime`, `Time2FaceFct`, `Time2Restore` |
| steady start | `rad_m1/rad_m1_clean.cpp` | `CleanSteadyStart` |
| restart state | `rad_m1/rad_m1_time2.cpp` | `Time2RstPack`, `Time2RstSet` |
| dt_retry checkpoint | `driver/dt_retry.cpp` | `DtRetry` |
| layout-invariant sums | `utils/repro_sum.hpp`, `rad_m1/rad_m1_layreduce.hpp`, `mesh/mesh.cpp` | `ReproSum`, `M1LayReduce` |

The files `rad_m1_clean_p1.cpp`, `_p2.cpp` and `rad_m1_clean_closure_p1.cpp`, `_p2.cpp` only
re-include their parent file in parts (`M1_TU_PART`) to shorten GPU compile times.

---

## 11. How to cite

If you use this scheme, please cite:

- **AthenaK**: J. M. Stone et al. (2024), *AthenaK: a performance-portable version of the
  Athena++ AMR framework*, ApJS. (TODO: insert the final volume and page numbers or DOI, and
  check the arXiv number, 2409.16053.)
- **This radiation scheme**: TODO, the method paper for the clean staggered M1–VET scheme
  (not yet published). Until it appears, cite this document with the commit hash of the code
  you ran.

Methods the scheme builds on, which should be cited where they are relevant:

- Ma, Pakmor, Justham & de Mink (2025), arXiv:2503.16627: the mixed-frame β terms (eqs. 56–57).
- Jiang (2021), ApJS 253, 49: discrete-ordinates implicit radiation in Athena++; the
  level-symmetric sets and the scattering split.
- Davis, Stone & Jiang (2012), ApJS 199, 9: VET closure from short characteristics in Athena.
- Krumholz, Klein, McKee & Bolstad (2007), ApJ 667, 626: mixed-frame source terms and the
  $O(\beta^2)$ terms in dynamic diffusion.
- Alexander (1977), SIAM J. Numer. Anal. 14, 1006: the L-stable two-stage SDIRK.
- Zalesak (1979), J. Comput. Phys. 31, 335: flux-corrected transport limiter.
- Eisenstat & Walker (1996), SIAM J. Sci. Comput. 17, 16: inexact Newton forcing.
- van der Vorst (1992), SIAM J. Sci. Stat. Comput. 13, 631: BiCGStab.

---

## 12. Limitations and open points

- **Uniform meshes only and $\hat c=c$** (`CleanInit` is FATAL otherwise).
- **Momentum on the streaming share** is not exactly conserved where `clean_stag_fwb` or
  `clean_stag_fcone` act. The imbalance is reported, not exchanged with the gas (§7.2). In
  3-D spherical polar with `clean_stag_fwb_ring = 0` (the production setting), the fwb residual's
  $L_z$ goes to the diagnostic ledger and is not projected out.
- **Layout invariance** is exact for the sums in every mode. The multigrid preconditioner is
  MeshBlock-local outside `tile` mode, so the converged answer can depend on the MeshBlock size
  at about the solver tolerance (`clean_conv_max`). In `tile` mode, the round-off of the
  cross-MeshBlock joins of the F line solves may remain layout-dependent. That statement comes
  from the development notes and was not re-measured for this document.
- **Code defaults versus production values.** Several keys are on by default for a fresh run
  but switched off in the reference production decks: `clean_stag_thinf`, `clean_stag_fwb_ring`,
  `clean_stag_lat_sums`, `clean_stag_opac_drag`, `clean_stag_tie_gid` and rotation. Others are
  off by default but on in production: `clean_stag_central`, `clean_stag_upwind_central`,
  `clean_face_hctau` and `clean_stag_eicarry`. Name these keys explicitly in production inputs.
- **Two separate face fluxes.** The face energy flux $\Phi_f$ that moves $E$ is the eliminated
  form of §4.2. The face momentum $F_f$ that the gas feels comes from the F solve of §4.6, which
  adds a thick-cell dissipation and the streaming-share transport of $F$. They agree where those
  extra terms vanish. Their difference is a discretisation choice: energy and momentum are each
  conserved in their own flux form. The precise relation between the two at a general fixed
  point is not documented in the code comments.
- **Inert keys.** `vet_source`, `vet_scatter_j` and the other relaxed-source keys in the
  production decks have no effect on the clean rays unless `clean_ray_source_relaxed = true`.
  `vet_col_lat_taucut` is overridden: clean rays sweep from the first shell, and the depth cut is
  `clean_ray_cut`'s.
