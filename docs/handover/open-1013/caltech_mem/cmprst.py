import sys
def data(f):
    b=open(f,'rb').read(); i=b.find(b'<par_end>'); return b[b.index(b'\n',i)+1:]
a,b=sys.argv[1],sys.argv[2]
da,db=data(a),data(b)
print('BITWISE' if da==db else 'DIFF (len %d vs %d)'%(len(da),len(db)), a, b)
