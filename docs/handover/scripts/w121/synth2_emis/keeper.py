"""List / download files from the pRT Keeper share (Seafile API; pRT's own fetcher needs Chrome)."""
import json
import os
import sys
import urllib.parse
import urllib.request
S = 'https://keeper.mpdl.mpg.de/api/v2.1/share-links/ccf25082fda448c8a0d0/dirents/?path='
D = 'https://keeper.mpdl.mpg.de/d/ccf25082fda448c8a0d0/files/?dl=1&p='
ROOT = '/viper/ptmp2/jinma/prt_data/input_data'


def ls(path):
    with urllib.request.urlopen(S + urllib.parse.quote(path)) as r:
        return json.load(r)['dirent_list']


def walk(path):
    for e in ls(path):
        if e['is_dir']:
            yield from walk(e['folder_path'])
        else:
            yield path.rstrip('/') + '/' + e['file_name'], e['size']


def get(path):
    out = ROOT + path
    if os.path.exists(out):
        return out
    os.makedirs(os.path.dirname(out), exist_ok=True)
    urllib.request.urlretrieve(D + urllib.parse.quote(path), out + '.part')
    os.rename(out + '.part', out)
    return out


if __name__ == '__main__':
    if sys.argv[1] == 'ls':
        for p, s in walk(sys.argv[2]):
            print('%10.1f MB  %s' % (s/1e6, p))
    else:
        for p in sys.argv[2:]:
            print(get(p), flush=True)
