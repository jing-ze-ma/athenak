#!/usr/bin/env python3
"""Compare the DATA PAYLOAD of two AthenaK .bin dumps -- everything after the embedded
parameter dump, whose LENGTH differs whenever the two runs carry different <problem>
keys (or different values of the same key).  The header states its own size."""
import sys


def payload(path):
    b = open(path, 'rb').read()
    i = b.find(b'  header offset=')
    j = b.find(b'\n', i)
    n = int(b[i + len(b'  header offset='):j])
    return b[j + 1 + n:]


a, b = payload(sys.argv[1]), payload(sys.argv[2])
print('IDENTICAL' if a == b else 'DIFFER (%d vs %d bytes)' % (len(a), len(b)))
