import sys
def pay(f):
    b=open(f,'rb').read(); i=b.index(b'<par_end>'); return b[i+len('<par_end>'):]
a=pay(sys.argv[1]); b=pay(sys.argv[2])
print(len(a),len(b), 'BITWISE' if a==b else 'DIFFER')
if a!=b:
    n=min(len(a),len(b)); d=[i for i in range(n) if a[i]!=b[i]]
    print('ndiff bytes',len(d),'first',d[:5])
