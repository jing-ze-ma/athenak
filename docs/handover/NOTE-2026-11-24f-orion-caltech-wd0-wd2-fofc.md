# NOTE 2026-11-24f orion -> Caltech: t9l62 WD0 vs WD2 (FOFC on), branch caltech-t9floor-wd2-140c271e0

**Binary.** Branch `caltech-t9floor-wd2-140c271e0` = vet-ia123-rc1 da92efe70 + three commits:
- 89a72f78f floor_audit_worst (device gather);
- 8197e0c4f implicit_pos_floor never floors a non-finite E;
- LAYC2 a2b44d59e (sdirk2: the u0 ghosts hold their neighbour's cell after Time2GasTime(3), so the flux limiters' donors agree on both blocks).

Build CUDA 12.9 + g++ 13 and run the tiny gate. Then run the 1-rank suite tests floor_pay + inf_reject on one H200, as in 24a.

**Arms.** Both from t9l62.00012.rst with `WD_t9l62_restart.args` (every key named), `hydro/floor_audit_worst=true` (already in the args), `time/nlim=6480`, `job/final_restart=false`, rst/bin dt 1e30. Same binary for both; same node or two nodes.
- **WD0**: the args as they are (FOFC off). This is the same-binary reference.
- **WD2**: the args + `hydro/fofc=true hydro/fofc_retest=true`. The IA11 rst stores fofc = false (USER 11-06, production setting). This 25-cycle diagnostic turns it on from the command line. Check that the log echo shows `fofc = true`.

**Report per arm** (to caltech_results/IA121_T9L62/arm_WD0/ and arm_WD2/: grep.txt + run.log + hst):
1. The `[update]` floor_audit line per cycle: floor_de, skip_de, and their per-cycle differences.
2. The `floor_audit_worst` lines in the blob (rank 0 gid 0, kji 76-80 48-49 47-57; r 2.29-2.45e12). In WD these created 1e41-2e42 per update at c6460-6478.
3. FOFC: the `FOFC first-order cell-stages` total (EOS floors line at the end) and the `fofc_retest: retested / unflagged` numbers in the floor_audit pay line.
4. CLEAN WARNING / stage fallback / DT_RETRY counts, and dt (min, last) over the 25 cycles.

**What to watch for (11-21 EY, WBFAST t9h):** FOFC on with the clean scheme gave CW 0 -> 9, hydro dt 350 -> 0.012 and floor_pay x18, while the log said `FOFC first-order cell-stages=0`. That points to an FOFC / clean interaction, or a counter that misses its own replacements.
- If WD2 shows a dt collapse or more CW while the first-order count reads 0, stop WD2 and send the log. It is the EY bug, not the answer.

**How we read it:**
- WD2's blob creation gone (worst dE << 1e40) and the first-order count > 0 in that blob: the high-order flux drives the deficit. FOFC is the existing by-construction fix; whether to use it in production is the USER's call.
- WD2 = WD0: the drain is in the hydro source terms (WB / geometric / am_transport), not the flux. Orion instruments those next.

-- orion (coordinator)
