#!/bin/bash
# resusage.sh <src file relative to src/> : device-only compile of the WORKTREE copy with
# kernel-resource-usage remarks; prints kernel (demangled, shortened), VGPRs, scratch
module purge >/dev/null 2>&1; module load gcc/14 rocm/6.3 openmpi_gpu/5.0 >/dev/null 2>&1
WT=/viper/ptmp2/jinma/wt_fast4; B=/viper/ptmp2/jinma/fast4_0925/b_m6_none_gpu; S=/viper/ptmp2/jinma/fast4_0925/src_m6
F=$1; O=/viper/ptmp2/jinma/fast4_0925/ru_$(basename $F .cpp)${TAGRU}.txt
hipcc -DKOKKOS_DEPENDENCE -DUSE_PROF_API=1 -D__HIP_PLATFORM_AMD__=1 -DMPI_PARALLEL_ENABLED=1 \
  -I$WT/src -I$WT -I$B -I$B/kokkos -I$B/kokkos/core/src -I$S/kokkos/core/src \
  -I$S/kokkos/tpls/desul/include -I$S/kokkos/tpls/mdspan/include -I$B/kokkos/containers/src \
  -I$S/kokkos/containers/src -I$B/kokkos/algorithms/src -I$S/kokkos/algorithms/src \
  -I$B/kokkos/simd/src -I$S/kokkos/simd/src \
  -isystem /viper/u2/system/soft/RHEL_9/packages/znver4/openmpi_gpu/gcc_14-14.1.0-rocm_6.3-6.3.4/5.0.8/include \
  -O3 -DNDEBUG -std=c++17 -fno-gpu-rdc --offload-arch=gfx942 -x hip --offload-device-only \
  -Rpass-analysis=kernel-resource-usage $XF -c $WT/src/$F -o /dev/null > $O.raw 2>&1
python3 - $O.raw <<'PY' > $O
import re,sys,subprocess
rows=[]; cur=None
for ln in open(sys.argv[1]):
    m=re.search(r"Function Name: (\S+)",ln)
    if m:
        cur=[m.group(1),{}]; rows.append(cur); continue
    m=re.search(r"remark:\s+(VGPRs|ScratchSize \[bytes/lane\]|Occupancy \[waves/SIMD\]|VGPRs Spill|Dynamic Stack): (\w+)",ln)
    if m and cur: cur[1][m.group(1).split()[0]+('S' if 'Spill' in m.group(1) else '')]=m.group(2)
dm=subprocess.run(['/usr/bin/c++filt'],input='\n'.join(r[0] for r in rows),capture_output=True,text=True).stdout.split('\n')
for (n,d),x in zip(rows,dm):
    m=re.search(r'RadiationM1::(\w+)\(',x); f=m.group(1) if m else re.sub(r'\(.*','',x)[-50:]
    l=re.findall(r'\{lambda\(([^)]*)\)#(\d+)\}',x)
    s=f+(''.join('#%s/%d'%(b,a.count('int')) for a,b in l[:1]))+(' R' if 'ParallelReduce' in x else '')
    if int(d.get('ScratchSize','0'))>0 or 'RadiationM1' in x or 'M1' in x:
        print('VGPR %4s spill %4s scratch %5s occ %2s dyn %5s  %s'%(d.get('VGPRs'),d.get('VGPRsS'),d.get('ScratchSize'),d.get('Occupancy'),d.get('Dynamic'),s))
PY
echo "done: $O ($(wc -l < $O) kernels with scratch or M1)"
