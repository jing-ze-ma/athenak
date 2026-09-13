# `lhlld` vs `hlld`: CPAW convergence + Orszag-Tang null test

CPU binary `build/src/athena` (serial, viper). All runs `rk2` + `plm`, ideal EOS, gamma = 5/3.

## A. Circularly polarized Alfven wave (3D oblique, travelling)

`inputs/tests/cpaw3d.athinput`, defaults (`b_par=1`, `b_perp=0.1`, `v_par=0`, `along_x1/2/3=false`),
`tlim=1.0`. Error = column 4 (0-based) of `CPAW_*-errs.dat`, the RMS of the per-variable L1 errors.

| solver | L1 (32x16x16) | L1 (64x32x32) | ratio L1_hi/L1_lo | order log2(L1_lo/L1_hi) | wall (lo/hi) |
| --- | --- | --- | --- | --- | --- |
| `hlld`  | 1.679649e-02 | 4.885990e-03 | 0.2909 | 1.78 | 0.49 s / 6.5 s |
| `lhlld` | 1.679297e-02 | 4.885848e-03 | 0.2909 | 1.78 | 0.53 s / 7.3 s |

Literature expectation: a CPAW is an exact nonlinear solution, so a 2nd-order scheme must show
2nd-order L1 convergence (ratio 0.25, order 2) for **both** solvers; the low-Mach (Minoshima &
Miyoshi) fix must not degrade the order. **PASS**: ratio 0.291 <= ~0.3 for both, and the two
solvers agree to ~2e-5 relative — the fix is essentially inactive here (this wave is not low-Mach:
|u| ~ b_perp/sqrt(rho) = 0.1 vs c_f ~ 1).

Command lines (run in separate scratch dirs; `RS` in {hlld,lhlld}, `(N,n2)` in {(32,16),(64,32)}):

```
athena -i inputs/tests/cpaw3d.athinput \
  job/basename=CPAW_${RS}_${N} \
  mesh/nx1=$N mesh/nx2=$n2 mesh/nx3=$n2 \
  meshblock/nx1=$N meshblock/nx2=$n2 meshblock/nx3=$n2 \
  time/tlim=1.0 time/integrator=rk2 mhd/reconstruct=plm mhd/rsolver=$RS \
  output1/dt=-1 output2/dt=-1 output3/dt=-1
```

## B. Orszag-Tang null test (t = 0.5)

`inputs/mhd/orszag_tang.athinput`, rk2 + plm, CFL 0.4, binary dump of `mhd_w_bcc` at t = 0.5 only.

| quantity | value |
| --- | --- |
| mean rho (128^2, hlld) | 0.221049 |
| d1 = <\|rho_lhlld,128 - rho_hlld,128\|> | 3.297137e-04 |
| d2 = <\|rho_hlld,128 - restrict(rho_hlld,256)\|> | 4.128856e-03 |
| **d1 / d2** | **0.0799** |
| fraction of cells with \|u\| < c_f/2  (phi < ~0.75) | 0.577 |
| fraction of cells with \|u\| < c_f/10 | 0.103 |
| max \|u\|/c_f in the domain | 1.86 |

`c_f = sqrt(0.5*((cs^2+ca^2) + sqrt((cs^2+ca^2)^2 - 4 cs^2 cax^2)))`, `cax` from `bcc1`,
`cs^2 = gamma (gamma-1) eint / rho`, all from the 128^2 `lhlld` dump at t = 0.5.
The 256^2 field is restricted to 128^2 by 2x2 block averaging.

Literature expectation: Minoshima & Miyoshi's supersonic/transonic tests (Orszag-Tang among them)
require the low-Mach correction to be a *null* change where the flow is transonic or supersonic —
the all-speed factor phi -> 1 there, so the scheme must reduce to plain HLLD. **PASS**:
d1/d2 = 0.08 << 1, i.e. the lhlld-vs-hlld difference is ~8% of hlld's own 128->256 discretization
error, well inside "order 0.1 or below". The phi proxy shows the fix is nevertheless touching a
real part of the domain (58% of cells below c_f/2, 10% below c_f/10), so the null result is not
simply because the correction never activates — the shocked/transonic regions where it is switched
off are the ones that dominate the solution error.

Command lines:

```
# 128^2, RS in {hlld,lhlld}
athena -i inputs/mhd/orszag_tang.athinput \
  job/basename=OT_${RS}_128 mesh/nx1=128 mesh/nx2=128 \
  meshblock/nx1=64 meshblock/nx2=64 time/tlim=0.5 mhd/rsolver=$RS \
  output1/dt=10.0 output2/file_type=bin output2/variable=mhd_w_bcc output2/dt=0.5 output3/dt=-1

# 256^2 reference
athena -i inputs/mhd/orszag_tang.athinput \
  job/basename=OT_hlld_256 mesh/nx1=256 mesh/nx2=256 \
  meshblock/nx1=128 meshblock/nx2=128 time/tlim=0.5 mhd/rsolver=hlld \
  output1/dt=10.0 output2/file_type=bin output2/variable=mhd_w_bcc output2/dt=0.5 output3/dt=-1
```

Dumps read with `vis/python/bin_convert.read_binary`, MeshBlocks reassembled onto the root grid
via `mb_logical` (i,j) and `nx1_out_mb`/`nx2_out_mb`.
