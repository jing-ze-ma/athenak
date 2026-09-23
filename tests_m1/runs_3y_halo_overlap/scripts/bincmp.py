# bincmp.py a b: 0 if the data after the parameter dump (<par_end>) are identical
import sys
def body(p):
    b = open(p, 'rb').read(); i = b.find(b'<par_end>')
    return b[i:] if i >= 0 else b
sys.exit(0 if body(sys.argv[1]) == body(sys.argv[2]) else 1)
