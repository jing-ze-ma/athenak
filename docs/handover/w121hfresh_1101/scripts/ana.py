# w121hfresh_1101 analysis: python3 ana.py on_a off_a off_b on_b
import re, sys, numpy as np
D = '/viper/ptmp2/jinma/w121hfresh_1101/runs/'
P = 110154.0
for r in sys.argv[1:]:
    log = open(D + r + '/run.log').read()
    rows = [tuple(map(float, m)) for m in
            re.findall(r'elapsed=(\S+) cycle=(\d+) time=(\S+) dt=(\S+)', log)]
    a = np.array(rows)
    # steady part: drop first 2 diag rows (startup + first 100 cycles)
    b = a[2:]
    spc = (b[-1, 0] - b[0, 0]) / (b[-1, 1] - b[0, 1])
    dtm = (b[-1, 2] - b[0, 2]) / (b[-1, 1] - b[0, 1])
    cpr = P / dtm
    nf = len(re.findall('FATAL', log))
    ev = [l.split() for l in open(D + r + '/dhj.log') if not l.startswith('#')]
    ev = np.array([[float(x) for x in e] for e in ev])
    mem = open(D + r + '/mem.log').read()
    tot = [int(x) for x in re.findall(r'MemTotal:\s+(\d+)', mem)]
    av = [int(x) for x in re.findall(r'MemAvailable:\s+(\d+)', mem)]
    used = (np.array(tot) - np.array(av)) / 1048576.0
    gtt = [int(x) for x in re.findall(r'GTT Total Used Memory \(B\):\s+(\d+)', mem)]
    print(f"{r}: cycles {int(a[-1,1])} t_end {a[-1,2]:.4e} s ({a[-1,2]/P:.4f} rot) "
          f"s/cycle {spc:.4f} mean dt {dtm:.3f} s last dt {a[-1,3]:.3f} cyc/rot {cpr:.0f} "
          f"wall/rot 1node {cpr*spc/60:.1f} min FATAL {nf}")
    print(f"   node mem used (MemTotal-MemAvail) start {used[0]:.1f} GiB peak {used.max():.1f} GiB"
          f" -> run {used.max()-used[0]:.1f} GiB; GTT peak per GPU {max(gtt)/2**30 if gtt else 0:.2f} GiB")
    print("   events (sum over log rows, cols dfloor efloor tfloor vceil fail c2p_it fofc):",
          ev[:, 1:8].sum(0).astype(int).tolist(), "rows", len(ev))
    h = np.loadtxt(D + r + '/dhj.user.hst')
    hd = open(D + r + '/dhj.user.hst').readlines()[1]
    lab = re.findall(r'\[\d+\]=(\S+)', hd)
    for c in ['Lrel', 'Iw', 'Tcor', 'Tcc', 'Tcont', 'Tbnd', 'Tcons', 'Twall', 'Efloor', 'Mfloor']:
        if c in lab:
            j = lab.index(c)
            print(f"   {c:7s} first {h[1,j]: .4e} mid {h[len(h)//2,j]: .4e} last {h[-1,j]: .4e}")
