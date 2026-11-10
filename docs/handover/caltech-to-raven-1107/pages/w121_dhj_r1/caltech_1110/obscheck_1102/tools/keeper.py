import json, sys, subprocess, os, urllib.parse
T='ccf25082fda448c8a0d0'; ROOT='/resnick/groups/carnegie_poc/jingze/obscheck_1102/prt_data/input_data'
def ls(path):
    out=subprocess.run(['curl','-s',f'https://keeper.mpdl.mpg.de/api/v2.1/share-links/{T}/dirents/?path='+urllib.parse.quote(path)],capture_output=True,text=True).stdout
    return [(e.get('folder_name') or e.get('file_name'), e.get('is_dir', 'folder_name' in e), e.get('size')) for e in json.loads(out)['dirent_list']]
def walk(path, depth=3):
    for n,isd,s in ls(path):
        p=path.rstrip('/')+'/'+n
        if isd and depth>0: yield from walk(p, depth-1)
        elif not isd: yield p, s
def get(p):
    dst=ROOT+p; os.makedirs(os.path.dirname(dst),exist_ok=True)
    if os.path.exists(dst): return dst
    subprocess.run(['curl','-sL','-o',dst+'.part',f'https://keeper.mpdl.mpg.de/d/{T}/files/?p='+urllib.parse.quote(p)+'&dl=1'],check=True)
    os.rename(dst+'.part',dst); return dst
if __name__=='__main__':
    if sys.argv[1]=='ls':
        for p,s in walk(sys.argv[2], int(sys.argv[3]) if len(sys.argv)>3 else 3): print(s,p)
    else:
        for p in sys.argv[2:]: print(get(p), flush=True)
