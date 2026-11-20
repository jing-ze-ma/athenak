#!/usr/bin/env python3
# ncusum.py DIR: one line per profiled kernel from k<i>.csv (ncu --page raw --csv): time, DRAM GB/s and % of peak, SM %, achieved/
# theoretical occupancy, occupancy limiter, regs, FP64 pipe %, L1/L2 hit, global load sectors per request (coalescing: 4 = ideal for
# 8-byte loads... see note), local-memory (spill) bytes, warp stall top reasons.
import sys, glob, csv, os, re
def val(row, *names):
    for n in names:
        if n in row and row[n] not in ('', 'n/a'):
            try: return float(row[n].replace(',', ''))
            except ValueError: return row[n]
    return float('nan')
top = {}
tf = os.path.join(sys.argv[1], 'top.txt')
if os.path.exists(tf):
    for i, l in enumerate([l for l in open(tf) if not l.startswith('#')]): top[f'k{i+1}'] = l.split('\t')[0]
print(f'{"k":4} {"key":34} {"us":>8} {"DRAM GB/s":>9} {"%pk":>5} {"SM%":>5} {"occA":>5} {"occT":>5} {"reg":>4} {"fp64%":>6} {"L1hit":>6} {"L2hit":>6} {"sec/req":>7} {"lmemB":>10}  stalls')
for f in sorted(glob.glob(os.path.join(sys.argv[1], 'k*.csv')), key=lambda x: int(re.findall(r'k(\d+)', x)[-1])):
    rows = list(csv.reader(l for l in open(f, errors='ignore') if l.startswith('"')))
    if len(rows) < 3: print(os.path.basename(f), 'no data'); continue
    h = rows[0]; r = dict(zip(h, rows[2])); U = dict(zip(h, rows[1]))
    SC = {'nsecond': 1e-3, 'usecond': 1.0, 'msecond': 1e3, 'second': 1e6, 'byte/second': 1e-9, 'Kbyte/second': 1e-6,
          'Mbyte/second': 1e-3, 'Gbyte/second': 1.0, 'Tbyte/second': 1e3, 'byte': 1.0, 'Kbyte': 1e3, 'Mbyte': 1e6, 'Gbyte': 1e9}
    for n in h:
        if U.get(n) in SC:
            try: r[n] = str(float(r[n].replace(',', ''))*SC[U[n]])
            except ValueError: pass
    k = os.path.basename(f)[:-4]
    t = val(r, 'gpu__time_duration.sum')
    dr = val(r, 'dram__bytes.sum.per_second', 'dram__bytes_read.sum.per_second')
    dp = val(r, 'dram__throughput.avg.pct_of_peak_sustained_elapsed', 'gpu__dram_throughput.avg.pct_of_peak_sustained_elapsed')
    sm = val(r, 'sm__throughput.avg.pct_of_peak_sustained_elapsed')
    oa = val(r, 'sm__warps_active.avg.pct_of_peak_sustained_active'); ot = val(r, 'sm__maximum_warps_per_active_cycle_pct')
    rg = val(r, 'launch__registers_per_thread')
    f64 = val(r, 'sm__pipe_fp64_cycles_active.avg.pct_of_peak_sustained_active', 'sm__inst_executed_pipe_fp64.avg.pct_of_peak_sustained_active')
    l1 = val(r, 'l1tex__t_sector_hit_rate.pct'); l2 = val(r, 'lts__t_sector_hit_rate.pct')
    sq = val(r, 'l1tex__t_sectors_pipe_lsu_mem_global_op_ld.sum'); rq = val(r, 'l1tex__t_requests_pipe_lsu_mem_global_op_ld.sum')
    lm = 0.0
    for n in ('l1tex__t_bytes_pipe_lsu_mem_local_op_ld.sum', 'l1tex__t_bytes_pipe_lsu_mem_local_op_st.sum'):
        v = val(r, n); lm += v if isinstance(v, float) and v == v else 0.0
    sn = [n for n in h if 'issue_stalled_' in n and n.endswith('per_issue_active.ratio')] or \
         [n for n in h if n.startswith('smsp__pcsamp_warps_issue_stalled_') and not n.endswith('not_issued')]
    st = sorted(((val(r, n), re.sub(r'.*issue_stalled_|_per_issue_active.ratio', '', n)) for n in sn), reverse=True)
    st = [(v, n) for v, n in st if isinstance(v, float) and v == v][:3]
    def g(x, fm): return fm % x if isinstance(x, float) else str(x)[:6]
    print(f'{k:4} {top.get(k, "?")[:34]:34} {g(t, "%8.0f")} {g(dr, "%9.0f")} {g(dp, "%5.1f")} {g(sm, "%5.1f")} {g(oa, "%5.1f")} {g(ot, "%5.1f")} {g(rg, "%4.0f")} {g(f64, "%6.1f")} {g(l1, "%6.1f")} {g(l2, "%6.1f")} {g(sq/rq if isinstance(sq, float) and isinstance(rq, float) and rq else float("nan"), "%7.1f")} {g(lm, "%10.0f")}  ' + ' '.join(f'{n}:{v:.1f}' for v, n in st))
