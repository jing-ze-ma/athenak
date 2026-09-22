"""Comparison 2: vertex-band vs panel-interior tangential flow, last dump of each
reconstruction arm (and of the restart arms, for reference).  Read-only."""
import sys
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
import bin_convert                       # noqa: E402
import dhjcs                             # noqa: E402
import common as C                       # noqa: E402

NB = 4        # vertex band half-width in cells


def panel_index(raw):
    """global (J,K) index inside the panel, 0..nx2_panel-1, for every (m,k,j)."""
    g = np.asarray(raw['mb_geometry'])
    nmb = g.shape[0]
    nj, nk = raw['nx2_mb'], raw['nx3_mb']
    n2 = raw['Nx2']
    n3 = raw['Nx3']
    J = np.zeros((nmb, nk, nj), dtype=int)
    K = np.zeros((nmb, nk, nj), dtype=int)
    for m in range(nmb):
        j0 = int(round((g[m, 2] - (-1.0)) / 2.0 * n2))
        k0 = int(round((g[m, 4] - (-1.0)) / 2.0 * n3))
        J[m] = (j0 + np.arange(nj))[None, :]
        K[m] = (k0 + np.arange(nk))[:, None]
    return J, K, n2, n3


def masks(raw):
    J, K, n2, n3 = panel_index(raw)
    nearJ = (J < NB) | (J >= n2 - NB)
    nearK = (K < NB) | (K >= n3 - NB)
    vert = nearJ & nearK
    inter = (J >= n2 // 4) & (J < 3 * n2 // 4) & (K >= n3 // 4) & (K < 3 * n3 // 4)
    edge = (nearJ ^ nearK)          # seam band but not a vertex
    return vert, inter, edge


def residual_speed(raw):
    """|v_tangential| and its departure from the zonal mean at the same latitude
    band and radius (the part a grid artefact lives in)."""
    geo = dhjcs.cs_geometry(raw)
    u, v = dhjcs.winds(raw, geo)        # zonal, meridional (physical, orthonormal)
    lat = geo['lat']
    vt = np.sqrt(u * u + v * v)
    lats, um = dhjcs.zonal_mean(np.moveaxis(u, 3, 0), lat)      # (ni, nlat)
    _, vm = dhjcs.zonal_mean(np.moveaxis(v, 3, 0), lat)
    il = np.clip(((np.degrees(lat) + 90) // 6).astype(int), 0, len(lats) - 1)
    du = u - np.moveaxis(um[:, il], 0, 3)
    dv = v - np.moveaxis(vm[:, il], 0, 3)
    return vt, np.sqrt(du * du + dv * dv), geo


def grid_noise(raw):
    """grid-scale (2 dx) roughness of the tangential velocity: the |1,-2,1| stencil in
    j and k, evaluated only where it fits INSIDE a MeshBlock (no halos in the dump)."""
    mb = raw['mb_data']
    v2 = np.asarray(mb['vely'], dtype=np.float64)
    v3 = np.asarray(mb['velz'], dtype=np.float64)
    n = np.full(v2.shape, np.nan)
    lj = np.abs(v2[:, :, 2:, :] - 2 * v2[:, :, 1:-1, :] + v2[:, :, :-2, :])
    lk = np.abs(v3[:, 2:, :, :] - 2 * v3[:, 1:-1, :, :] + v3[:, :-2, :, :])
    n[:, 1:-1, 1:-1, :] = np.sqrt(lj[:, 1:-1, :, :]**2 + lk[:, :, 1:-1, :]**2)
    return n


def one(path):
    raw = bin_convert.read_binary(path)
    vert, inter, edge = masks(raw)
    vt, res, geo = residual_speed(raw)
    gn = grid_noise(raw)
    g = C.cs_metric(raw)
    out = dict(t=raw['time'], rot=raw['time'] / C.PROT, r=g['rc'] / 9.44e9)
    for nm, msk in [('vert', vert), ('inter', inter), ('edge', edge)]:
        out['vt_' + nm] = np.sqrt(np.mean(vt[msk]**2, axis=0))
        out['res_' + nm] = np.sqrt(np.mean(res[msk]**2, axis=0))
        out['gn_' + nm] = np.sqrt(np.nanmean(gn[msk]**2, axis=0))
        out['n_' + nm] = int(msk.sum())
    return out


def main():
    fh = open(C.OUT + 'vertex_tables.md', 'w')
    fh.write('# Vertex band vs panel interior (vertex.py)\n\n'
             'Vertex band = within %d cells of a cube vertex in BOTH gnomonic angles; '
             'interior = the middle half of each panel in both. `res` = rms departure '
             'of the horizontal wind from its zonal mean at the same latitude and '
             'radius.\n' % NB)
    res = {}
    for group, arms in [('comparison 2 (cold start)', C.ARMS2),
                        ('comparison 1 (restart, reference)', C.ARMS1)]:
        fh.write('\n### %s\n\n| arm | rot | cells vert/inter | rms |v_t| vert | '
                 'interior | ratio | rms residual vert | interior | ratio | '
                 'grid-noise vert | interior | ratio |\n'
                 % group + '|---|---|---|---|---|---|---|---|---|---|---|---|\n')
        for name, d in arms:
            fs = sorted(glob.glob(d + '/bin/*.bin'))
            # comparison 2: the three arms stopped at different times, so take the dump
            # nearest the latest time they all reached (6.8e4 s = 0.22 rot).
            f = fs[-1] if arms is C.ARMS1 else min(
                fs[1:], key=lambda p: abs(bin_convert.read_binary(p)['time'] - 6.8e4))
            o = one(f)
            res[name] = o
            # mass-blind global rms over all radii
            a = np.sqrt(np.mean(o['vt_vert']**2))
            b = np.sqrt(np.mean(o['vt_inter']**2))
            c = np.sqrt(np.mean(o['res_vert']**2))
            e = np.sqrt(np.mean(o['res_inter']**2))
            p = np.sqrt(np.nanmean(o['gn_vert']**2))
            q = np.sqrt(np.nanmean(o['gn_inter']**2))
            fh.write('| %s | %.3f | %d/%d | %.3e | %.3e | %.3f | %.3e | %.3e | %.3f '
                     '| %.3e | %.3e | %.3f |\n'
                     % (name, o['rot'], o['n_vert'], o['n_inter'], a, b, a / b, c, e,
                        c / e, p, q, p / q))
    # radial detail for comparison 2
    fh.write('\n### comparison 2: vertex/interior residual ratio vs radius\n\n'
             '| r/Rp | ' + ' | '.join(n for n, _ in C.ARMS2) + ' |\n|---|'
             + '---|' * len(C.ARMS2) + '\n')
    r = res['plm_ctl']['r']
    for i in range(0, len(r), 12):
        fh.write('| %.3f | ' % r[i] + ' | '.join(
            '%.2f' % (res[n]['res_vert'][i] / res[n]['res_inter'][i])
            for n, _ in C.ARMS2) + ' |\n')
    fh.close()

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for n, _ in C.ARMS2:
        o = res[n]
        ax[0].semilogy(o['r'], o['res_vert'], '-', label=n + ' vertex')
        ax[0].semilogy(o['r'], o['res_inter'], '--', alpha=.6, label=n + ' interior')
        ax[1].plot(o['r'], o['res_vert'] / o['res_inter'], label=n)
    ax[0].set_ylabel('rms |v - zonal mean| [cm/s]')
    ax[1].set_ylabel('vertex / interior')
    ax[1].axhline(1, color='k', lw=.8)
    for a_ in ax:
        a_.set_xlabel('r / Rp')
        a_.grid(alpha=.3)
        a_.legend(fontsize=7)
    fig.suptitle('Comparison 2 (cold start, nghost=3): vertex-localised residual flow')
    fig.tight_layout()
    fig.savefig(C.OUT + 'fig_vertex.png', dpi=130)
    print(open(C.OUT + 'vertex_tables.md').read())


main()
