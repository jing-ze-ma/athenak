import re
p = '/viper/ptmp2/jinma/wt_accretor/src/coordinates/coordinates.cpp'
s = open(p).read()
pat = re.compile(r"(?P<lead>[ \t]*)(?P<decl>(?:const )?Real pr = )gen_ \? wder_\(m,IDPR,(?P<ix>[a-z, ]+)\)\n"
                 r"[ \t]*: eos_\.Pressure\(w0\(m,IDN,(?P=ix)\), w0\(m,IEN,(?P=ix)\)\);")
def rep(mo):
    lead, decl, ix = mo.group('lead'), mo.group('decl'), mo.group('ix')
    pad = ' ' * (len(lead) + len(decl))
    return (f"{lead}// isothermal EOS: no energy in w0 (IEN is out of range), p = d c_s^2\n"
            f"{lead}{decl}(!eos_.is_ideal) ? w0(m,IDN,{ix})*SQR(eos_.iso_cs)\n"
            f"{pad}: gen_ ? wder_(m,IDPR,{ix})\n"
            f"{pad}: eos_.Pressure(w0(m,IDN,{ix}), w0(m,IEN,{ix}));")
s2, n = pat.subn(rep, s)
print('replaced', n)
open(p, 'w').write(s2)
