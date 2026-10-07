import pypdf,glob,os,sys
os.chdir(sys.argv[1])
for f in glob.glob('*.pdf'):
    t=f[:-4]+'.txt'
    if os.path.exists(t): continue
    try:
        r=pypdf.PdfReader(f); open(t,'w').write('\n'.join((p.extract_text() or '') for p in r.pages))
    except Exception as e: print(f,e)
