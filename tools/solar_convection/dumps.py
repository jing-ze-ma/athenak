"""Locate a run's binary dumps and reject the truncated ones.

A GPFS per-user quota can silently truncate an AthenaK .bin dump mid-write: the run
carries on, the file is short, and nothing in the log says so. Reading one gives
garbage (bin_convert happily returns whatever meshblocks made it to disk, the rest
uninitialised), which shows up as nonsense frames in the movies.

Every dump of a given run has the same size to within a couple of bytes (only the ASCII
header numbers vary), so a file more than TOL bytes below the run's largest dump is
truncated. `dump_files` returns a list with the SAME indexing as the on-disk sequence,
with truncated entries replaced by None so frame numbers keep matching file numbers.
"""

import glob
import os

TOL = 65536      # bytes a healthy dump may fall below the largest one


def dump_files(run_dir, pattern='sun.hydro_w.*.bin', verbose=True):
    """Sorted dump paths for one run, truncated ones replaced by None."""
    files = sorted(glob.glob(os.path.join(run_dir, 'bin', pattern)))
    if not files:
        raise RuntimeError('no dumps matching %s in %s' % (pattern, run_dir))
    sizes = [os.path.getsize(f) for f in files]
    full = max(sizes)
    out = [f if s >= full - TOL else None for f, s in zip(files, sizes)]
    bad = [(i, f, s) for i, (f, s) in enumerate(zip(files, sizes)) if out[i] is None]
    if bad and verbose:
        print('WARNING: %d of %d dumps in %s are truncated (expected %d bytes) and '
              'will be skipped:' % (len(bad), len(files), run_dir, full))
        for i, f, s in bad:
            print('    [%d] %s  %d bytes (%.1f%%)'
                  % (i, os.path.basename(f), s, 100.0*s/full))
    return out


def last_good(files):
    """Index of the last non-truncated dump."""
    for i in range(len(files) - 1, -1, -1):
        if files[i] is not None:
            return i
    raise RuntimeError('every dump is truncated')


def resolve(files, idx):
    """Turn an index into the index of an intact dump.

    Non-negative indices address the on-disk sequence directly (so they keep matching
    the file numbers) and must not be truncated. Negative indices count back through
    the INTACT dumps only, so -1 is always the last usable snapshot.
    """
    good = [i for i, f in enumerate(files) if f is not None]
    if not good:
        raise RuntimeError('every dump is truncated')
    if idx < 0:
        if -idx > len(good):
            raise IndexError('snapshot %d out of range (%d intact dumps)'
                             % (idx, len(good)))
        return good[idx]
    if idx >= len(files):
        raise IndexError('snapshot %d out of range (%d dumps)' % (idx, len(files)))
    if files[idx] is None:
        raise RuntimeError('snapshot %d is truncated' % idx)
    return idx
