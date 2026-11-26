"""Plot style of the AREPO-IDORT Zenodo notebooks (record 15032138): mesa2.mplstyle + black text.
   import idort_style; idort_style.use()          # line plots: font 16 (wave_2d: font=20)
"""
import os
import matplotlib as mpl
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


def use(font=None):
    mpl.style.use('default')
    plt.style.use(os.path.join(HERE, 'mesa2.mplstyle'))
    for k in ('text.color', 'axes.labelcolor', 'axes.edgecolor', 'xtick.color', 'ytick.color'):
        mpl.rcParams[k] = 'black'
    mpl.rcParams['figure.facecolor'] = 'white'
    if font:
        mpl.rcParams['font.size'] = font
