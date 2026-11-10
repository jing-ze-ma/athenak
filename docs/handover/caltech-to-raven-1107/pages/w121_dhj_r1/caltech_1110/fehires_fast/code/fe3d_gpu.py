"""fe3d_gpu.py: GPU (CuPy) path of fe3d.py, same physics and same outputs.

Per exposure the ray geometry, the KD-tree lookup and the per-sample interpolation run on the CPU
(vectorised, NPROC worker processes, one exposure each, ahead of the GPU); the GPU does, per batch
of rays and for every pressure cut at once:
  - the E_low x velocity line-profile matrices Gm, Am, Bm (one batched DGEMM; samples outside the
    per-cut mask lw > max - 25 get weight 0, i.e. are left out exactly as on the CPU),
  - the per-line peak test (tau_peak > taumin),
  - the line deposit + excess absorption, as a per-pixel GATHER: each pixel thread adds the lines
    in the CPU order (line index ascending; per line the q term then the 1-q term), compiled
    without FMA contraction, then exp(-tau_c)(1 - exp(-tau)) * area is added to the spectrum row
    in the CPU ray order.  Differences to the CPU path: libm vs CUDA exp/expm1 ulps, the DGEMM
    summation order, the order of the per-ray template-weighted sum W.
Called by fe3d.run() when gpu=1 (or env FEHIRES_GPU=1); env FEHIRES_GPU_BATCH (rays per batch, 256).
"""
import os
import time
import numpy as np
import fe3d as M
from fe3d import G, P, F

BATCH = int(os.environ.get('FEHIRES_GPU_BATCH', '256'))
MAXC = 8
GATHER2 = os.environ.get('FEHIRES_GATHER', '1') == '2'   # 11-03: no-atomics W (per-block partials)

KSRC = r'''
// NC (number of cuts) is a compile-time constant (-DNC=..): the accumulators stay in registers
__device__ __forceinline__ double tline(const double* __restrict__ Gb, const long long rs,
    const int nj, const int nE, const int e, const double f, const double cc, const double a,
    const double bb, const int j) {
  const double* g0 = Gb + (long long)e*rs + j;
  const double* g1 = g0 + rs;
  const double* b0 = Gb + (long long)(nE + e)*rs + nj + j;
  const double* b1 = b0 + rs;
  return cc*((g0[0]*(1 - f) + g1[0]*f) + a*(g0[nj]*(1 - f) + g1[nj]*f)
             + bb*(b0[0]*(1 - f) + b1[0]*f));
}
extern "C" __global__ void fe_gather(
    const int B, const int nc, const int nl, const int npix, const int nv, const int nE,
    const long long* __restrict__ iq, const double* __restrict__ fq, const int* __restrict__ ie,
    const double* __restrict__ fe, const double* __restrict__ lc, const double* __restrict__ gr,
    const double* __restrict__ gw, const int* __restrict__ klo, const int* __restrict__ khi,
    const unsigned char* __restrict__ keep, const double* __restrict__ Out,
    const double* __restrict__ ltauc, const int* __restrict__ nodei,
    const double* __restrict__ nodew, const double* __restrict__ tmpl,
    const double* __restrict__ area, const int* __restrict__ limb,
    const unsigned char* __restrict__ act, double* __restrict__ spec, double* __restrict__ W)
{
  const int p = blockIdx.x*blockDim.x + threadIdx.x;
  const bool in = p < npix;
  const int nj = 2*nv + 1;
  const long long rs = 2*nj;
  const long long bs = (long long)2*nc*nE*rs;
  // warp: the union of the lines that reach any of its 32 pixels (lines sorted by pixel); every
  // lane evaluates the profile at its own j once; the q term of pixel p is lane p-1's value
  const int lane = threadIdx.x & 31;
  const int pw = p - lane;
  const int k0 = pw < npix ? klo[pw] : 0;
  const int k1 = pw < npix ? khi[min(pw + 31, npix - 1)] : 0;
  double acc0[NC], acc1[NC];
  int ni = 0; double nw = 0.0, tm = 0.0;
#pragma unroll
  for (int c = 0; c < NC; c++) { acc0[c] = 0.0; acc1[c] = 0.0; }
  if (in) {
    ni = nodei[p]; nw = nodew[p]; tm = tmpl[p];
#pragma unroll
    for (int c = 0; c < NC; c++) {
      acc0[c] = spec[(2*c)*(long long)npix + p];
      acc1[c] = spec[(2*c + 1)*(long long)npix + p];
    }
  }
  for (int b = 0; b < B; b++) {
    const int lb = limb[b];
#pragma unroll
    for (int c = 0; c < NC; c++) {
      if (!act[b*nc + c]) continue;               // uniform over the grid
      double ex = 0.0;
      {                                              // all lanes: warp shuffles inside
        const double* Gb = Out + b*bs + (long long)(2*c)*nE*rs;
        const unsigned char* kp = keep + ((long long)b*nc + c)*nl;
        double tau = 0.0;
        for (int k = k0; k < k1; k++) {
          if (!kp[k]) continue;                      // uniform over the warp
          const int j = p - (int)(iq[k] - nv);
          const int e = ie[k];
          const double f = fe[k], cc = lc[k], a = gr[k], bb = gw[k], q = fq[k];
          const bool vj = in && j >= 0 && j < nj;
          const double tj = vj ? tline(Gb, rs, nj, nE, e, f, cc, a, bb, j) : 0.0;
          double tq = __shfl_up_sync(0xffffffff, tj, 1);
          const bool vm = in && j - 1 >= 0 && j - 1 < nj;
          if (lane == 0) tq = vm ? tline(Gb, rs, nj, nE, e, f, cc, a, bb, j - 1) : 0.0;
          if (vm && tq > 1e-14) tau += tq*q;
          if (vj && tj > 1e-14) tau += tj*(1 - q);
        }
        if (in && tau != 0.0 && tau > 1e-10) {
          const double* lt_ = ltauc + b*9;
          const double lt = lt_[ni]*(1 - nw) + lt_[ni + 1]*nw;
          ex = exp(-exp(lt))*(-expm1(-tau))*area[b];
          if (lb) acc1[c] += ex; else acc0[c] += ex;
        }
      }
      double wv = ex*tm;
      if (__any_sync(0xffffffff, wv != 0.0)) {
        for (int o = 16; o > 0; o >>= 1) wv += __shfl_down_sync(0xffffffff, wv, o);
        if ((threadIdx.x & 31) == 0 && wv != 0.0) atomicAdd(W + b*nc + c, wv);
      }
    }
  }
  if (in) {
#pragma unroll
    for (int c = 0; c < NC; c++) {
      spec[(2*c)*(long long)npix + p] = acc0[c];
      spec[(2*c + 1)*(long long)npix + p] = acc1[c];
    }
  }
}

extern "C" __global__ void fe_gather2(
    const int B, const int nc, const int nl, const int npix, const int nv, const int nE,
    const long long* __restrict__ iq, const double* __restrict__ fq, const int* __restrict__ ie,
    const double* __restrict__ fe, const double* __restrict__ lc, const double* __restrict__ gr,
    const double* __restrict__ gw, const int* __restrict__ klo, const int* __restrict__ khi,
    const unsigned char* __restrict__ keep, const double* __restrict__ Out,
    const double* __restrict__ ltauc, const int* __restrict__ nodei,
    const double* __restrict__ nodew, const double* __restrict__ tmpl,
    const double* __restrict__ area, const int* __restrict__ limb,
    const unsigned char* __restrict__ act, double* __restrict__ spec, double* __restrict__ W)
{
  extern __shared__ double wsh[];                // [nwarp][B*nc]
  const int wid = threadIdx.x >> 5, nwarp = blockDim.x >> 5;
  for (int q = threadIdx.x; q < nwarp*B*nc; q += blockDim.x) wsh[q] = 0.0;
  __syncthreads();
  const int p = blockIdx.x*blockDim.x + threadIdx.x;
  const bool in = p < npix;
  const int nj = 2*nv + 1;
  const long long rs = 2*nj;
  const long long bs = (long long)2*nc*nE*rs;
  // warp: the union of the lines that reach any of its 32 pixels (lines sorted by pixel); every
  // lane evaluates the profile at its own j once; the q term of pixel p is lane p-1's value
  const int lane = threadIdx.x & 31;
  const int pw = p - lane;
  const int k0 = pw < npix ? klo[pw] : 0;
  const int k1 = pw < npix ? khi[min(pw + 31, npix - 1)] : 0;
  double acc0[NC], acc1[NC];
  int ni = 0; double nw = 0.0, tm = 0.0;
#pragma unroll
  for (int c = 0; c < NC; c++) { acc0[c] = 0.0; acc1[c] = 0.0; }
  if (in) {
    ni = nodei[p]; nw = nodew[p]; tm = tmpl[p];
#pragma unroll
    for (int c = 0; c < NC; c++) {
      acc0[c] = spec[(2*c)*(long long)npix + p];
      acc1[c] = spec[(2*c + 1)*(long long)npix + p];
    }
  }
  for (int b = 0; b < B; b++) {
    const int lb = limb[b];
#pragma unroll
    for (int c = 0; c < NC; c++) {
      if (!act[b*nc + c]) continue;               // uniform over the grid
      double ex = 0.0;
      {                                              // all lanes: warp shuffles inside
        const double* Gb = Out + b*bs + (long long)(2*c)*nE*rs;
        const unsigned char* kp = keep + ((long long)b*nc + c)*nl;
        double tau = 0.0;
        for (int k = k0; k < k1; k++) {
          if (!kp[k]) continue;                      // uniform over the warp
          const int j = p - (int)(iq[k] - nv);
          const int e = ie[k];
          const double f = fe[k], cc = lc[k], a = gr[k], bb = gw[k], q = fq[k];
          const bool vj = in && j >= 0 && j < nj;
          const double tj = vj ? tline(Gb, rs, nj, nE, e, f, cc, a, bb, j) : 0.0;
          double tq = __shfl_up_sync(0xffffffff, tj, 1);
          const bool vm = in && j - 1 >= 0 && j - 1 < nj;
          if (lane == 0) tq = vm ? tline(Gb, rs, nj, nE, e, f, cc, a, bb, j - 1) : 0.0;
          if (vm && tq > 1e-14) tau += tq*q;
          if (vj && tj > 1e-14) tau += tj*(1 - q);
        }
        if (in && tau != 0.0 && tau > 1e-10) {
          const double* lt_ = ltauc + b*9;
          const double lt = lt_[ni]*(1 - nw) + lt_[ni + 1]*nw;
          ex = exp(-exp(lt))*(-expm1(-tau))*area[b];
          if (lb) acc1[c] += ex; else acc0[c] += ex;
        }
      }
      double wv = ex*tm;
      if (__any_sync(0xffffffff, wv != 0.0)) {
        for (int o = 16; o > 0; o >>= 1) wv += __shfl_down_sync(0xffffffff, wv, o);
        if ((threadIdx.x & 31) == 0) wsh[wid*B*nc + b*nc + c] = wv;
      }
    }
  }
  if (in) {
#pragma unroll
    for (int c = 0; c < NC; c++) {
      spec[(2*c)*(long long)npix + p] = acc0[c];
      spec[(2*c + 1)*(long long)npix + p] = acc1[c];
    }
  }
  __syncthreads();
  for (int q = threadIdx.x; q < B*nc; q += blockDim.x) {   // W here = per-block partials [block][B*nc]
    double t = 0.0;
    for (int w = 0; w < nwarp; w++) t += wsh[w*B*nc + q];
    W[(long long)blockIdx.x*B*nc + q] = t;
  }
}
'''


def _geom(ph, kv):
    """the ray list of one exposure, in the CPU order (theta outer, b inner)."""
    a = 2*np.pi*ph
    ob = np.array([np.cos(a), -np.sin(a), 0.0])
    e1 = np.array([np.sin(a), np.cos(a), 0.0])
    e2 = np.array([0.0, 0.0, 1.0])
    dth = float(kv.get('dth', 2.0))
    th = np.deg2rad(np.arange(0, 360, dth) + dth/2)
    ds = float(kv.get('ds', 5e7))
    bb, db = G['bb'], G['db']
    r1 = G['r1']
    xs, ys = M.AOR*np.sin(a), M.BIMP*np.cos(a)
    u1, u2 = P['limb_darkening']
    ldn = 1 - u1/3 - u2/6
    R = dict(it=[], ib=[], limb=[], area=[], ds=[], r=[], pts=[])
    for it, t_ in enumerate(th):
        limb = 0 if np.cos(t_) > 0 else 1
        sx = xs + bb*np.cos(t_)/M.RS
        sy = ys + bb*np.sin(t_)/M.RS
        onstar = (sx**2 + sy**2) < 1.0
        if u1 != 0.0 or u2 != 0.0:
            mu1 = 1 - np.sqrt(np.clip(1 - sx**2 - sy**2, 0, 1))
            ldw = (1 - u1*mu1 - u2*mu1**2)/ldn
        else:
            ldw = None
        for ib, b_ in enumerate(bb):
            if not onstar[ib]:
                continue
            smax = np.sqrt(max(r1**2 - b_**2, 0.0))
            ns = max(int(2*smax/ds), 4)
            s = (np.arange(ns) + 0.5)/ns*2*smax - smax
            P0 = b_*(np.cos(t_)*e1 + np.sin(t_)*e2)
            pts = P0[None, :] + s[:, None]*ob[None, :]
            r = np.linalg.norm(pts, axis=1)
            ok = (r >= G['rc'][0]) & (r <= r1)
            if ok.sum() < 2:
                continue
            area = b_*db[ib]*np.deg2rad(dth)
            if ldw is not None:
                area = area*ldw[ib]
            R['it'].append(it)
            R['ib'].append(ib)
            R['limb'].append(limb)
            R['area'].append(area)
            R['ds'].append(2*smax/ns)
            R['r'].append(r[ok])
            R['pts'].append(pts[ok]/r[ok][:, None])
    return R, ob, th


SHM = None                                  # /dev/shm dir for worker -> GPU process arrays


def prep_phase(args):
    """CPU part of one exposure (worker process): rays and the KD-tree column lookup; the
    per-sample arrays go through SHM (np.save), the interpolation is done on the GPU."""
    ph, kv = args
    t0 = time.time()
    R, ob, th = _geom(ph, kv)
    out = dict(ph=ph, ob=ob, th=th, it=np.array(R['it'], int), ib=np.array(R['ib'], int),
               limb=np.array(R['limb'], np.int32), area=np.array(R['area'], float),
               ds=np.array(R['ds'], float), n=np.array([len(x) for x in R['r']], int))
    out['offs'] = np.concatenate([[0], np.cumsum(out['n'])]).astype(int)
    if R['r']:
        r = np.concatenate(R['r'])
        _, col = G['tree'].query(np.concatenate(R['pts']))
        rc = G['rc']
        fi = np.interp(r, rc, np.arange(len(rc)))
        i0 = np.clip(fi.astype(int), 0, len(rc) - 2)
        f = fi - i0
        fn = os.path.join(SHM, 'ph%d_%.6f.npz' % (os.getpid(), ph))
        np.savez(fn, col=col.astype(np.int32), i0=i0.astype(np.int32), f=f)
        out['file'] = fn
    out['sec'] = time.time() - t0
    return out


FUSEK = r'''
extern "C" __global__ void interp(const long long Sb, const int S, const int nc, const int nr,
    const long long ncell, const int* __restrict__ col, const int* __restrict__ i0,
    const double* __restrict__ f_, const long long* __restrict__ rr, const long long* __restrict__ jj,
    const double* __restrict__ dsr, const double ob0, const double ob1,
    const float* __restrict__ AT, const float* __restrict__ AVx, const float* __restrict__ AVy,
    const float* __restrict__ AnP, const float* __restrict__ Alp, const float* __restrict__ Alw,
    const float* __restrict__ Alac, const int B,
    double* Tp, double* vp, double* nPp, double* lpp, double* lwp, double* ex9) {
  // ip(A) = A[c, i0]*(1 - f) + A[c, i0 + 1]*f, as fe3d.ray_one; outputs scattered to (B, S) pads
  const long long i = (long long)blockIdx.x*blockDim.x + threadIdx.x;
  if (i >= Sb) return;
  const long long a0 = (long long)col[i]*nr + i0[i], a1 = a0 + 1;
  const double f = f_[i], g = 1 - f;
  const long long o = rr[i]*S + jj[i];
  const double T = (double)AT[a0]*g + (double)AT[a1]*f;
  Tp[o] = T;
  vp[o] = -(((double)AVx[a0]*g + (double)AVx[a1]*f)*ob0 + ((double)AVy[a0]*g + (double)AVy[a1]*f)*ob1);
  nPp[o] = ((double)AnP[a0]*g + (double)AnP[a1]*f)*pow(T/1e4, 0.3);
  lpp[o] = (double)Alp[a0]*g + (double)Alp[a1]*f;
  for (int c = 0; c < nc; c++)
    lwp[(long long)c*B*S + o] = (double)Alw[c*ncell + a0]*g + (double)Alw[c*ncell + a1]*f;
  const double ds = dsr[rr[i]];
  for (int n = 0; n < 9; n++)
    ex9[o*9 + n] = exp((double)Alac[a0*9 + n]*g + (double)Alac[a1*9 + n]*f)*ds;
}
extern "C" __global__ void profR(const long long N, const int nj, const double* __restrict__ Tp,
    const double* __restrict__ vp, const double* __restrict__ VG, const double c2kb,
    const double mfe, const double sqpi, double* Rm) {
  // Rm[b, s, j]: Phi (j < nj) and Kw (j >= nj), as fe3d.ray_one
  const long long i = (long long)blockIdx.x*blockDim.x + threadIdx.x;
  if (i >= N) return;
  const long long bs = i / (2*nj);
  const int j = i % (2*nj);
  const double T = Tp[bs];
  const double bth = sqrt(c2kb*T/mfe);
  const double dv = VG[j < nj ? j : j - nj] - vp[bs];
  const double x = dv/bth;
  if (j < nj) Rm[i] = exp(-x*x)/(sqpi*bth);
  else Rm[i] = (1 - exp(-pow(x/2.5, 4.0)))/fmax(dv*dv, 1.0);
}
extern "C" __global__ void buildL(const int B, const int S, const int nE, const int nc,
    const double* __restrict__ Tp, const double* __restrict__ nPp, const double* __restrict__ w,
    const double* __restrict__ EG, const double ev, const double kb, double* L) {
  // L[b, 2c nE + e, s] = Bz, L[b, (2c+1) nE + e, s] = Bz nP;  Bz = exp(-(EV/(KB T)) E_e) w_c
  const long long i = (long long)blockIdx.x*blockDim.x + threadIdx.x;
  if (i >= (long long)B*nE*S) return;
  const int s = i % S;
  const int e = (i / S) % nE;
  const int b = i / ((long long)S*nE);
  const long long o = (long long)b*S + s;
  const double ez = exp(-((ev/(kb*Tp[o]))*EG[e]));
  const double np_ = nPp[o];
  for (int c = 0; c < nc; c++) {
    const double bz = ez*w[(long long)c*B*S + o];
    L[((long long)b*2*nc*nE + (2*c)*nE + e)*S + s] = bz;
    L[((long long)b*2*nc*nE + (2*c + 1)*nE + e)*S + s] = bz*np_;
  }
}
'''


SEGK = r'''
extern "C" __global__ void segsum(const int B, const int S, const int K, const long long* n,
                                  const double* x, double* out) {
  // out[b, k] = sum_{s < n_b} x[b, s, k], sequential in s (the CPU order)
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i >= B*K) return;
  const int b = i / K, k = i % K;
  double acc = 0.0;
  for (long long s = 0; s < n[b]; s++) acc += x[((long long)b*S + s)*K + k];
  out[i] = acc;
}
'''


class GPU:
    def __init__(self):
        import cupy as cp
        self.cp = cp
        nv, npix = M.NV, M.NPIX
        nj = 2*nv + 1
        iq = G['iq'].astype(np.int64)
        nl = len(iq)
        # pixel p receives from line k iff 0 <= p - (iq_k - nv) <= 2 nv + 1; lines sorted by pixel
        assert np.all(np.diff(iq) >= 0)
        lgl = np.log(F.LAM0) + F.DLN*np.arange(npix)
        lgn = np.log(M.LNODE)
        nodei = np.clip(np.searchsorted(lgn, lgl) - 1, 0, len(lgn) - 2)
        nodew = (lgl - lgn[nodei])/(lgn[nodei + 1] - lgn[nodei])
        a = cp.asarray
        self.d = dict(iq=a(iq), fq=a(G['fq']), ie=a(G['ie'].astype(np.int32)), fe=a(G['fe']),
                      lc=a(G['lc']), gr=a(G['gr']), gw=a(G['gw']),
                      klo=a(np.searchsorted(iq, np.arange(npix) - nv - 1, 'left').astype(np.int32)),
                      khi=a(np.searchsorted(iq, np.arange(npix) + nv, 'right').astype(np.int32)),
                      nodei=a(nodei.astype(np.int32)),
                      nodew=a(nodew), tmpl=a(np.ascontiguousarray(G['tmpl'], float)),
                      VG=a(M.VG), EG=a(M.EG), iel=a(G['ie'].astype(np.int64)))
        for k in ('lp', 'lac', 'T', 'Vx', 'Vy', 'nP', 'lw'):
            self.d[k] = a(G[k])
        self.nl = nl
        self.prof = {} if os.environ.get('FEHIRES_PROF') else None
        self.seg = cp.RawKernel(SEGK, 'segsum', options=('--fmad=false',))
        self.kint = cp.RawKernel(FUSEK, 'interp', options=('--fmad=false',))
        self.kR = cp.RawKernel(FUSEK, 'profR', options=('--fmad=false',))
        self.kL = cp.RawKernel(FUSEK, 'buildL', options=('--fmad=false',))
        self._kg = {}

    def kgather(self, nc):
        if nc not in self._kg:
            self._kg[nc] = self.cp.RawKernel(KSRC, 'fe_gather', options=('--fmad=false', '-DNC=%d' % nc))
        return self._kg[nc]

    def kgather2(self, nc):
        if ('g2', nc) not in self._kg:
            self._kg[('g2', nc)] = self.cp.RawKernel(KSRC, 'fe_gather2', options=('--fmad=false', '-DNC=%d' % nc))
        return self._kg[('g2', nc)]

    def tick(self, name):
        if self.prof is None:
            return
        self.cp.cuda.Device().synchronize()
        t = time.time()
        self.prof[name] = self.prof.get(name, 0.0) + t - self._t
        self._t = t

    def phase(self, pr, kv):
        cp = self.cp
        if self.prof is not None:
            cp.cuda.Device().synchronize()
            self._t = time.time()
        d = self.d
        t0 = time.time()
        taumin = float(kv.get('taumin', M.TAUL_MIN))
        nv, npix, nE = M.NV, M.NPIX, len(M.EG)
        nj = 2*nv + 1
        nc = G['lw'].shape[0]
        assert nc <= MAXC
        nray = len(pr['n'])
        spec = cp.zeros((nc, 2, npix))
        Wd = cp.zeros((nray, nc))
        sq = np.sqrt(np.pi)
        if nray and 'col_d' in pr:                       # prep=gpu: arrays already on the GPU
            col_a, i0_a, f_a = pr['col_d'], pr['i0_d'], pr['f_d']
        elif nray:
            z = np.load(pr['file'])
            col_a, i0_a, f_a = cp.asarray(z['col']), cp.asarray(z['i0']), cp.asarray(z['f'])
            os.remove(pr['file'])
        ob = pr['ob']
        lptv = np.zeros(nray)
        for b0 in range(0, nray, BATCH):
            b1 = min(b0 + BATCH, nray)
            B = b1 - b0
            nn = pr['n'][b0:b1]
            S = int(nn.max())
            s0, s1 = pr['offs'][b0], pr['offs'][b1]
            rr_h = np.repeat(np.arange(B), nn)
            rr = cp.asarray(rr_h)
            jj = cp.asarray(np.arange(s1 - s0) - (pr['offs'][b0:b1] - s0)[rr_h])
            c_, i_, f = col_a[s0:s1], i0_a[s0:s1], f_a[s0:s1]
            dsr = cp.asarray(pr['ds'][b0:b1])

            Sb = int(s1 - s0)
            Tp = cp.full((B, S), 1000.0)
            vp = cp.zeros((B, S))
            nPp = cp.zeros((B, S))
            lpp = cp.full((B, S), -np.inf)
            lwp = cp.full((nc, B, S), -np.inf)
            ex9 = cp.zeros((B, S, 9))
            ncell = d['T'].size
            self.kint(((Sb + 127)//128,), (128,),
                      (np.int64(Sb), np.int32(S), np.int32(nc), np.int32(d['T'].shape[1]),
                       np.int64(ncell), c_, i_, f, rr, jj, dsr, np.float64(ob[0]),
                       np.float64(ob[1]), d['T'], d['Vx'], d['Vy'], d['nP'], d['lp'], d['lw'],
                       d['lac'], np.int32(B), Tp, vp, nPp, lpp, lwp, ex9))
            lpt = lpp.max(1)
            tauc = cp.empty((B, 9))
            nn_d = cp.asarray(nn.astype(np.int64))
            self.seg(((B*9 + 127)//128,), (128,), (np.int32(B), np.int32(S), np.int32(9), nn_d,
                                                   ex9, tauc))
            del ex9
            tok = tauc.min(1) <= 15
            lwm = lwp.max(2)                                               # (nc, B)
            acc_ = tok[None, :] & (lwm >= -600)
            act = cp.ascontiguousarray(acc_.T).astype(cp.uint8)            # (B, nc)
            w = cp.where((lwp > (lwm - 25)[:, :, None]) & acc_[:, :, None],
                         cp.exp(lwp)*dsr[None, :, None], 0.0)
            del lwp
            ltauc = cp.log(cp.maximum(tauc, 1e-300))
            lptv[b0:b1] = lpt.get()
            self.tick('interp')
            Rm = cp.empty((B, S, 2*nj))
            NR = B*S*2*nj
            self.kR(((NR + 127)//128,), (128,),
                    (np.int64(NR), np.int32(nj), Tp, vp, d['VG'], np.float64(2*F.KB),
                     np.float64(F.M_FE), np.float64(sq), Rm))
            L = cp.empty((B, 2*nc*nE, S))
            NL_ = B*nE*S
            self.kL(((NL_ + 127)//128,), (128,),
                    (np.int32(B), np.int32(S), np.int32(nE), np.int32(nc), Tp, nPp, w,
                     d['EG'], np.float64(F.EV), np.float64(F.KB), L))
            self.tick('profiles')
            Out = cp.ascontiguousarray(cp.matmul(L, Rm))                   # (B, 2 nc nE, 2 nj)
            self.tick('dgemm')
            del L, Rm
            keep = cp.empty((B, nc, self.nl), cp.uint8)
            for c in range(nc):
                gmax = Out[:, 2*c*nE:(2*c + 1)*nE, :nj].max(2)              # (B, nE)
                tpk = d['lc'][None, :]*(gmax[:, d['iel']]*(1 - d['fe'][None, :])
                                        + gmax[:, d['iel'] + 1]*d['fe'][None, :])
                keep[:, c, :] = (tpk > taumin) & (act[:, c:c + 1] > 0)
            if self.prof is not None:
                self.prof['nkeep'] = self.prof.get('nkeep', 0.0) + float(keep.sum())
                self.prof['nbc'] = self.prof.get('nbc', 0.0) + float(act.sum())
            self.tick('keep')
            if GATHER2 and B*nc*8*8 <= 48*1024:
                nt = 256
                nbk = (npix + nt - 1)//nt
                Wp = cp.empty((nbk, B, nc))
                self.kgather2(nc)((nbk,), (nt,),
                   (np.int32(B), np.int32(nc), np.int32(self.nl), np.int32(npix), np.int32(nv),
                    np.int32(nE), d['iq'], d['fq'], d['ie'], d['fe'], d['lc'], d['gr'], d['gw'],
                    d['klo'], d['khi'], keep, Out, ltauc, d['nodei'],
                    d['nodew'], d['tmpl'], cp.asarray(pr['area'][b0:b1]),
                    cp.asarray(pr['limb'][b0:b1]), act, spec, Wp), shared_mem=(nt//32)*B*nc*8)
                Wd[b0:b1] = Wp.sum(0)
                self.tick('gather')
                del Out, keep
                continue
            Wb = cp.zeros((B, nc))
            nt = 128
            self.kgather(nc)(((npix + nt - 1)//nt,), (nt,),
                   (np.int32(B), np.int32(nc), np.int32(self.nl), np.int32(npix), np.int32(nv),
                    np.int32(nE), d['iq'], d['fq'], d['ie'], d['fe'], d['lc'], d['gr'], d['gw'],
                    d['klo'], d['khi'], keep, Out, ltauc, d['nodei'],
                    d['nodew'], d['tmpl'], cp.asarray(pr['area'][b0:b1]),
                    cp.asarray(pr['limb'][b0:b1]), act, spec, Wb))
            Wd[b0:b1] = Wb
            self.tick('gather')
            del Out, keep
        th = pr['th']
        nb = len(G['bb'])
        Wm = np.zeros((nc, len(th), nb))
        Pt = np.zeros((len(th), nb))
        Wh = Wd.get()
        Wm[:, pr['it'], pr['ib']] = Wh.T
        Pt[pr['it'], pr['ib']] = lptv
        norm = np.pi*M.RS**2
        sp = spec.get()
        return dict(spec=(sp/norm).astype(np.float32), W=Wm/norm, ptan=Pt, th=th,
                    sec=pr['sec'] + time.time() - t0)


def run_phases(phases, kv):
    """all exposures: CPU prep in a Pool (forked before CUDA starts), GPU in this process."""
    from multiprocessing import Pool
    import shutil
    import tempfile
    global SHM
    SHM = tempfile.mkdtemp(prefix='fehgpu_', dir='/dev/shm' if os.path.isdir('/dev/shm') else None)
    out = []
    try:
        with Pool(min(M.NPROC, len(phases))) as pool:
            it = pool.imap(prep_phase, [(ph, kv) for ph in phases], chunksize=1)
            gpu = None
            for pr in it:
                if gpu is None:
                    gpu = GPU()
                t0 = time.time()
                out.append(gpu.phase(pr, kv))
                print('  phase %+.4f: prep %.1f s, gpu %.1f s, %d rays' % (
                    pr['ph'], pr['sec'], time.time() - t0, len(pr['n'])), flush=True)
    finally:
        shutil.rmtree(SHM, ignore_errors=True)
    return out


# ---------------------------------------------------------------------------------------------
# prep=gpu (11-03): the per-exposure ray geometry, the nearest-column lookup and the radial
# interpolation weights on the GPU, in the GPU process (no worker Pool, no /dev/shm files).
# Same arithmetic as _geom/prep_phase (separate elementwise kernels, --fmad=false, scalar cos/sin
# of theta on the CPU as before); the KD-tree is replaced by an exact brute-force nearest column
# (squared distance ((dx dx + dy dy) + dz dz), lowest index on ties): identical except exact ties.
PREPK = r'''
extern "C" __global__ void nncol(const long long N, const double* __restrict__ ux,
    const double* __restrict__ uy, const double* __restrict__ uz, const int nc,
    const double* __restrict__ XC, int* __restrict__ col) {
  extern __shared__ double sh[];
  const long long i = (long long)blockIdx.x*blockDim.x + threadIdx.x;
  double x = 0, y = 0, z = 0;
  if (i < N) { x = ux[i]; y = uy[i]; z = uz[i]; }
  double best = 1e300; int bi = 0;
  const int TILE = 1024;
  for (int t0 = 0; t0 < nc; t0 += TILE) {
    const int nt = min(TILE, nc - t0);
    __syncthreads();
    for (int k = threadIdx.x; k < 3*nt; k += blockDim.x) sh[k] = XC[3*t0 + k];
    __syncthreads();
    if (i < N) {
      for (int k = 0; k < nt; k++) {
        const double dx = x - sh[3*k], dy = y - sh[3*k + 1], dz = z - sh[3*k + 2];
        const double d = (dx*dx + dy*dy) + dz*dz;
        if (d < best) { best = d; bi = t0 + k; }
      }
    }
  }
  if (i < N) col[i] = bi;
}
extern "C" __global__ void interp_r(const long long N, const double* __restrict__ r,
    const int n, const double* __restrict__ xp, int* __restrict__ i0o, double* __restrict__ fo) {
  // np.interp(r, xp, arange(n)) as numpy's arr_interp (precomputed slopes), then
  // i0 = clip(int(fi), 0, n - 2), f = fi - i0
  const long long i = (long long)blockIdx.x*blockDim.x + threadIdx.x;
  if (i >= N) return;
  const double x = r[i];
  double fi;
  if (x < xp[0]) fi = 0.0;
  else if (x > xp[n - 1]) fi = (double)(n - 1);
  else if (x == xp[n - 1]) fi = (double)(n - 1);
  else {
    int lo = 0, hi = n - 1;                       // xp[lo] <= x < xp[hi]
    while (hi - lo > 1) { const int m = (lo + hi) >> 1; if (xp[m] <= x) lo = m; else hi = m; }
    if (xp[lo] == x) fi = (double)lo;
    else {
      const double slope = ((double)(lo + 1) - (double)lo)/(xp[lo + 1] - xp[lo]);
      fi = slope*(x - xp[lo]) + (double)lo;
    }
  }
  long long k = (long long)fi;
  k = k < 0 ? 0 : (k > n - 2 ? n - 2 : k);
  i0o[i] = (int)k;
  fo[i] = fi - (double)k;
}
'''


def gpu_prep(gpu, ph, kv):
    cp = gpu.cp
    t0 = time.time()
    if not hasattr(gpu, 'knn'):
        gpu.knn = cp.RawKernel(PREPK, 'nncol', options=('--fmad=false',))
        gpu.kir = cp.RawKernel(PREPK, 'interp_r', options=('--fmad=false',))
        gpu.XC = cp.asarray(np.ascontiguousarray(G['tree'].data, float).ravel())
        gpu.rcd = cp.asarray(np.ascontiguousarray(G['rc'], float))
    a = 2*np.pi*ph
    ob = np.array([np.cos(a), -np.sin(a), 0.0])
    e1 = np.array([np.sin(a), np.cos(a), 0.0])
    e2 = np.array([0.0, 0.0, 1.0])
    dth = float(kv.get('dth', 2.0))
    th = np.deg2rad(np.arange(0, 360, dth) + dth/2)
    ds = float(kv.get('ds', 5e7))
    bb, db = G['bb'], G['db']
    r1 = G['r1']
    xs, ys = M.AOR*np.sin(a), M.BIMP*np.cos(a)
    u1, u2 = P['limb_darkening']
    ldn = 1 - u1/3 - u2/6
    ct = np.array([np.cos(t_) for t_ in th])           # scalar calls, as the per-ray loop
    st = np.array([np.sin(t_) for t_ in th])
    sx = xs + bb[None, :]*ct[:, None]/M.RS
    sy = ys + bb[None, :]*st[:, None]/M.RS
    onstar = (sx**2 + sy**2) < 1.0
    nb = len(bb)
    area = np.empty((len(th), nb))
    for ib in range(nb):
        area[:, ib] = bb[ib]*db[ib]*np.deg2rad(dth)
    if u1 != 0.0 or u2 != 0.0:
        mu1 = 1 - np.sqrt(np.clip(1 - sx**2 - sy**2, 0, 1))
        ldw = (1 - u1*mu1 - u2*mu1**2)/ldn
        area = area*ldw
    smax = np.array([np.sqrt(max(r1**2 - b_**2, 0.0)) for b_ in bb])
    ns = np.array([max(int(2*sm/ds), 4) for sm in smax])
    SM = int(ns.max())
    sp = np.zeros((nb, SM))
    for ib in range(nb):
        sp[ib, :ns[ib]] = (np.arange(ns[ib]) + 0.5)/ns[ib]*2*smax[ib] - smax[ib]
    dsr = 2*smax/ns
    # P0 = b (cos t e1 + sin t e2) per (theta, b), component-wise as the 3-vector expression
    P0 = [bb[None, :]*(ct[:, None]*e1[k] + st[:, None]*e2[k]) for k in range(3)]
    spd = cp.asarray(sp)[None, :, :]
    X = cp.asarray(P0[0])[:, :, None] + spd*ob[0]
    Y = cp.asarray(P0[1])[:, :, None] + spd*ob[1]
    Z = cp.asarray(P0[2])[:, :, None] + spd*ob[2]
    r = cp.sqrt((X*X + Y*Y) + Z*Z)
    jv = cp.arange(SM)[None, None, :] < cp.asarray(ns)[None, :, None]
    ok = (r >= G['rc'][0]) & (r <= r1) & jv & cp.asarray(onstar)[:, :, None]
    cnt = ok.sum(2)
    valid = cnt >= 2
    ok &= valid[:, :, None]
    idx = cp.flatnonzero(ok)
    rs = r.ravel()[idx]
    ux = X.ravel()[idx]/rs
    uy = Y.ravel()[idx]/rs
    uz = Z.ravel()[idx]/rs
    del X, Y, Z, r, ok, jv
    N = int(idx.size)
    col = cp.empty(N, cp.int32)
    i0 = cp.empty(N, cp.int32)
    f = cp.empty(N, cp.float64)
    if N:
        gpu.knn(((N + 127)//128,), (128,), (np.int64(N), ux, uy, uz, np.int32(G['tree'].n), gpu.XC, col),
                shared_mem=3*1024*8)
        gpu.kir(((N + 127)//128,), (128,), (np.int64(N), rs, np.int32(len(G['rc'])), gpu.rcd, i0, f))
    vh = valid.get()
    it_, ib_ = np.nonzero(vh)                           # row-major = theta outer, b inner
    n = cnt.get()[it_, ib_].astype(int)
    pr = dict(ph=ph, ob=ob, th=th, it=it_.astype(int), ib=ib_.astype(int),
              limb=np.where(ct[it_] > 0, 0, 1).astype(np.int32), area=area[it_, ib_].astype(float),
              ds=dsr[ib_].astype(float), n=n, col_d=col, i0_d=i0, f_d=f)
    pr['offs'] = np.concatenate([[0], np.cumsum(n)]).astype(int)
    gpu.cp.cuda.Device().synchronize()
    pr['sec'] = time.time() - t0
    return pr


def run_phases_gpuprep(phases, kv):
    out = []
    gpu = GPU()
    for ph in phases:
        pr = gpu_prep(gpu, ph, kv)
        t0 = time.time()
        out.append(gpu.phase(pr, kv))
        print('  phase %+.4f: prep(gpu) %.2f s, gpu %.2f s, %d rays' % (
            pr['ph'], pr['sec'], time.time() - t0, len(pr['n'])), flush=True)
    if gpu.prof is not None:
        print('  prof', {k: round(v, 2) for k, v in gpu.prof.items()}, flush=True)
    return out
