"""shared figure style (v10 IDORT page: white surface, dark ink, recessive grid)"""
import numpy as np
import matplotlib.pyplot as plt
C1, C2, C3, C4 = '#2a78d6', '#eb6834', '#1baf7a', '#6250d6'
INK, INK2, GRID = '#1a1a1a', '#52514e', '#e4e4e2'
plt.rcParams.update({'font.size': 14, 'axes.titlesize': 14, 'axes.labelsize': 14, 'legend.fontsize': 11.5, 'xtick.labelsize': 13, 'ytick.labelsize': 13,
    'text.color': INK, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2, 'axes.edgecolor': '#9a9a96',
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.7, 'lines.linewidth': 1.8, 'figure.facecolor': 'white',
    'axes.facecolor': 'white', 'savefig.facecolor': 'white', 'legend.frameon': True, 'legend.framealpha': 0.92})
def slope2(ax, anchor=None):
    """reference slope -2 line through the first point of the first plotted line (or anchor (x,y))"""
    ls = [l for l in ax.get_lines() if len(l.get_xdata()) > 1]
    if not ls and anchor is None: return
    if anchor is None:
        x, y = np.asarray(ls[0].get_xdata(), float), np.asarray(ls[0].get_ydata(), float)
        x0, y0 = x[0], y[0] * 1.6; x1 = x[-1]
    else:
        (x0, y0), x1 = anchor[:2], anchor[2]
    xs = np.array([x0, x1]); ax.loglog(xs, y0 * (xs / x0) ** -2, color=INK, ls='--', lw=1.1, label='slope 2')
