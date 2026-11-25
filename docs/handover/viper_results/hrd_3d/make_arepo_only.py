# AREPO points only (no other codes, no other legend names), dark, normal background
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_hrd as m
rows = m.load(os.path.join(m.HERE, 'sims_3d.csv'))
m.draw(rows, True, 'hrd_3d_dark_arepo_only', dark=True, codes=['AREPO'], highlight=['AREPO'])
m.draw(rows, True, 'hrd_3d_dark_arepo_only_transparent', dark=True, transparent=True,
       codes=['AREPO'], highlight=['AREPO'])
