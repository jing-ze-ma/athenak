# w121hfresh_1101 1x: WASP-121b hydro fresh-start deck at 1x solar metallicity

`deck/hfresh_1x.athinput` is the 1x twin of `deck/hfresh.athinput` (3x, the tested fresh-start deck,
fork/fin51-r-fix 182477553 defaults). Every changed key carries `# hfresh_1x:` and its 3x value.
Nothing else changes: `diff` against hfresh.athinput shows only these lines and a 2-line header. The roche keys stay off,
and so does the tide.

## Keys changed (1x vs 3x)
| key | 1x | 3x | source |
|---|---|---|---|
| hydro/rad_met, problem/met | 0.0 | 0.4771 | 1x |
| hydro/eos_xh, eos_yhe | 0.7381, 0.2485 | 0.7188, 0.2420 | 1x solar; the same pair as the 10-01 1x deck (w121fresh_1x_f6). **Not in the brief's list, but needed:** the 3x EOS composition would make a 3x-EOS / 1x-opacity hybrid |
| problem/ck_data_dir | ck/CK1X | ck/CK3X | package-relative |
| problem/ck_table | ck/CK1X/ck/Premixed_1x_g8_11_hiT2.txt | ck/CK3X/ck/Premixed_3x_... | 1x table |
| problem/ck_ce_table | CE_tables/FastChem_ck_1x_int_hiT2.txt | ..._3x_... | 1x CE table |
| problem/ic_profile | ic/ic_w121_1x_f5.txt | ic/ic_w121_3x_f5.txt | 1x face-5 RCE IC (w121fresh_1001/ic, md5 10ee6d18...) |
| problem/ap = mesh/x1min | 1.122218e10 | 1.124210e10 | 1x anchor, see below |
| problem/grav | 1176.9609 | 1172.7924 | 1x anchor |
| mesh/x1max | 2.150381e10 | 2.1e10 | grid112 |
| mesh/f_stretch_r_c1..c8, b1/b2 amp/x/w | see deck | see deck | grid112 |
| nx1 (mesh, meshblock) | 112 (unchanged) | 112 | |

## Tables (`pkg/ck/CK1X/`)
- `ck/Premixed_1x_g8_11_hiT2.txt` md5 2df0d027b8f42e937091b66402a311a1, `CE_tables/FastChem_ck_1x_int_hiT2.txt` md5
  b6643409062f9ef3b070667c67b62da8. These are the repo files `data/exo_fms_ck/...` (untracked, mtime 09-25), the exact paths
  the 10-01 1x deck `w121fresh_1001/w121fresh_1x_f6.athinput` reads (its ck_data_dir is that repo directory).
  `ckimpl_budget_1027/data/CK1X` has no ck/ or CE_tables/ subdirectory, so it was not used.
- `cia/`, `ray/`, `sw_flux/` are identical to CK3X (md5 checked, all 12 files) and are **not** duplicated here.
  **To install: copy or symlink `pkg/ck/CK3X/{cia,ray,sw_flux}` into `pkg/ck/CK1X/`.** The deck reads `ck_data_dir/sw_flux/sw_band_flux_W121_11.txt`
  and the cia/ray files from ck_data_dir. Run directories link `ck` and `ic` to the package, as `job.sh` does
  (`job_1x.sh` here).

## Anchor: p_ref and ap / grav (`gate/`)
- p_ref (1x, NIRSpec NRS1 2.7-3.7 um) = **7.07e-4 bar** (wasp121_0925/x1/transit_1x.txt; 3x 3.956e-4, 10x 1.8e-4).
- Deck (tide off): planet_setup.py anchor_limb (the planet-package method, toi1518b/scripts copy = 3e0bd80a) on
  ic_w121_1x_f5.txt + the 1x EOS dump (deepconv_0925/dump/eos_table.txt, md5 697ca675...), M_p 1.170 M_J, R_p 1.742 R_J:
  **ap = x1min = 1.122218e10, grav = 1176.961** (`gate/planet_setup_1x_f5.txt`). roche_anchor.py --no_roche gives the same
  numbers to all digits. (The script calls np.trapezoid; the login-node numpy lacks it. It ran through a 3-line shim
  np.trapezoid = np.trapz with no other change.)
- The 10-01 1x deck's 1.126575e10 / 1167.8749 came from the older 1x IC ic_w121_1x.txt (wasp121_0925/x1). On the f5 IC the anchor
  is 0.39 % (4.36e7 cm) lower and grav 0.78 % higher. The 3x hfresh deck's 1.124210e10 / 1172.7924 was likewise derived on
  ic_w121_3x.txt (met3/setup_w121_3x.txt), not on the f5 IC that deck uses.
- With the tide (Caltech NOTE 272395a5b; roche_anchor.py from fork/fin51-r-fix c8aac483c), WASP-121b 1x:
  **ap = x1min = 1.144050e10, grav = 1132.4682, roche_a = 3.846161e11** (`gate/roche_anchor_1x.txt`, command line inside).
  The tool takes G M_* from omega and roche_a (Kepler), which implies M_* = 1.3938 M_sun vs the literature 1.358 (+2.63 %).
  - 3x regression (`gate/roche_anchor_3x_regress.txt`) with Caltech's command line (p_ref 1.8e-4, 3x EOS dump
    w121prod_0929/met3): ap 1.129091e10, grav 1162.6751, roche_a 3.846161e11. Caltech has 1.129094e10 / 1162.669, which
    agrees to 3e-6 / 5e-6, probably because they used a different 3x dump (eos_3x.txt).
  - **Caltech's 3x reference uses p_ref 1.8e-4 bar, which is the 10x NRS1 value.** At the 3x transit pressure 3.956e-4 bar
    the tool gives ap 1.141384e10, grav 1137.7664 (`gate/roche_anchor_3x_pref3x.txt`).
- A 1x grid at the Roche anchor (x1min = 1.144050e10) is in `grid/gridkeys_112_roche.txt` (same grid112 method:
  x1max 2.232001e10, max dr/target 1.058). It was NOT run or gated; with the tide x1max/r_L1 = 0.92.

## Radial grid (`grid/`, toi1518b/scripts/grid112.py, PLANET_EOS = 1x dump)
Cells per scale height of the 3x hfresh grid carried at equal pressure onto the 1x column. x1max = r(p(2.1e10 on the 3x
column) = 6.1e-18 bar) on the 1x column = 2.150381e10.
- **nx1 112 does not strictly meet the 3x cells/H everywhere:** the worst cell is 5.6 % coarser than the 3x n_H at its pressure
  (max dr/target 1.0558). Band minima, 1x vs 3x: 1e-6..1e-4 bar 3.09 vs 2.98, 1e-4..1e-2 3.95 vs 4.04, 1e-2..1 4.08 vs 4.10,
  1..10 3.75 vs 3.80, 10..300 3.96 vs 4.01. The min over p >= 1e-6 bar is 3.09 (3x 2.98). The target cell count is 112.1.
  nx1 108 gives a worst cell 9.8 % coarser. nx1 stays 112 (one MeshBlock in x1).
- Radial dt estimate (CFL 0.3): 11.72 s vs 11.43 s for the 3x (ratio 1.025).

## Gates
1. Parse (`-n`, CPU athena_dhj_cpu_8bcd50c8 = 182477553 + the remap merge; also GPU binary, job 12528605): rc 0. The keys
   are as listed. Unread-key warnings are the same 9 as the 3x deck.
2. CPU fresh start, 2 cycles, 12 ranks (`gate/gate.sh`, `gate.out`): rc 0, no FATAL, dt 14.317 s. The events at cycles 0/1/2
   are dfloor 1216/3648/3648 and efloor 472/1616/1804, all at the top (min density 1e-18; the 3x GPU smoke had 1744/984 at cycle 0).
   - IC T vs table (`gate/gate_ic_active.txt`, active cells with x1, x2 and x3 ghosts excluded, rho >= 1e-14, 581824 cells):
     max |dT/T| 3.9e-3 (1e-9..1e-6 bar), 8.5e-4 (1e-6..1e-3), 8.3e-5 (1e-3..1), 3.2e-5 (1..300), 1.5e-5 (300..3000).
     `gate_ic.txt` (toi1518b gate_ic.py, 1x EOS grid) reports 0.835 in 1e-6..1e-3 bar. Those are cubed-sphere panel-corner
     x2/x3 ghost cells, which that script counts as active. No active cell is off.
   - The deepest active cells lie on the table: 222.6 bar / 4230 K and 694.0 bar / 5006 K. The deepest inner ghost is 997 bar,
     rel -4.4e-6.
   - r(p_ref = 7.07e-4 bar) over 9600 columns is 1.232109e10 .. 1.259062e10, against the tool's pole 1.232081e10 and equator
     1.259071e10: **+2.3e-5 / -7.3e-6**.
3. GPU smoke (apudev, 1 node x 2 MI300A, athena_dhj_gpu72_18247755_hfresh md5 38c3fa6c..., HSA_XNACK=1,
   HSA_NO_SCRATCH_RECLAIM=1, `job_1x.sh`):
   - 40-cycle smoke, job 12528605: rc 0, no FATAL.
   - 10-min run, job 12528606 (vipa1001), from `ana.py` with the first 100 cycles excluded (`gate/ana_on1x_a.txt`):
     rc 0, **no FATAL**, 21700 cycles, 2.573 rot. Throughput **0.0272 s/cycle**, against 0.0277-0.0279 for the 3x ON arm
     (on_a/b/c, RESULTS.md).
     - **Mean dt 13.05 s** (3x: 13.69 s). dt does not settle as it did at 3x (13.79 s): it wanders between 12.15 and 13.5 s,
       12.15 s at the last row. That gives 8441 cycles per rotation (3x: 8044), **3.8 min per rotation** on 1 node (3x: 3.7).
     - Events summed over the dhj.log rows (1x vs 3x on_a): eos_dfloor 3.86e8 vs 4.87e8, eos_efloor 5.37e8 vs 1.10e9,
       eos_tfloor 2.10e7 vs 3.65e8. eos_vceil, eos_fail, c2p_it and fofc are 0 in both.
     - Memory: GTT peak 20.5 GiB per GPU and 59.7 GiB per node, the same as 3x.
4. md5s: `MD5SUMS`.
