"""tcmp.py run... : shell-mean v_r at the dumps of short tests (sponge off from rst 5), dv/dt by shell"""
import sys
src=open('/resnick/groups/carnegie_poc/jingze/r136noise_1106/an1.py').read().split('res={}')[0]
runs=sys.argv[1:];sys.argv=[sys.argv[0]]
exec(src.replace("RUN=R+'/runB/prod'","RUN='/resnick/groups/carnegie_poc/jingze/r136open_1106/runs/T0'"))
RR=(30,35,40,45,50,55,60,65,70,75,80,85,90,95,99.2)
idx=[np.argmin(abs(r-x*RSUN)) for x in RR]
print('shell r       '+' '.join('%7.1f'%(r[i]/RSUN) for i in idx))
print('dr            '+' '.join('%7.3f'%(dr[i]/RSUN) for i in idx))
OUT={}
for run in runs:
    RUN='/resnick/groups/carnegie_poc/jingze/r136open_1106/runs/'+run
    for n in (21,22):
        try: dh=bc.read_binary(RUN+'/bin/r136.hydro_w.%05d.bin'%n)
        except Exception as e: continue
        t=dh['time']-432000.;w=assemble(dh);v=w['velx'].astype(np.float64);rho=w['dens'].astype(np.float64)
        OUT.setdefault(run,{})[int(round(t))]=[float(v[...,i].mean()/t) for i in range(len(r))]
        print('%-4s t+%5.0f  <vr>/t '%(run,t)+' '.join('%7.0f'%(v[...,i].mean()/t) for i in idx)+'  rho-wt '+' '.join('%6.0f'%((rho*v)[...,i].sum()/rho[...,i].sum()/t) for i in idx[-5:]))
import json;json.dump(dict(r=list(map(float,r/RSUN)),dr=list(map(float,dr/RSUN)),acc=OUT),open('/resnick/groups/carnegie_poc/jingze/r136open_1106/tcmp.json','w'))
