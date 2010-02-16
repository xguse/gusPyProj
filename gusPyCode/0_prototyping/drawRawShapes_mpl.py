import matplotlib
#matplotlib.use('Qt4Agg')
matplotlib.use('Qt4Agg')
from matplotlib import pylab as pl


fig = pl.figure()
ax = fig.add_subplot(111)

zInc = 0

lines = []

lines.append(matplotlib.lines.Line2D((-2000,10),(1,1),color='black', linewidth=2, zorder=1))
lines.append(matplotlib.lines.Line2D((-2000,10),(2,2),color='black',linewidth=2, zorder=4))

for l in lines:
    ax.add_line(l)

patches = []

patches.append(matplotlib.patches.Rectangle( (3,0.75), label='MpRty', color='#EAFC71', width=20, height=0.5,  zorder=4))
patches.append(matplotlib.patches.Rectangle( (-1,1.75), width=3, height=0.5,  zorder=5))
patches.append(matplotlib.patches.Rectangle( (-789,1.75), label='MpRty', color='#EAFC71', width=20, height=0.5,  zorder=6))

for p in patches:
    ax.add_patch(p)



ax.autoscale_view()
#ax.figure.canvas.draw()
ax.set_xlim((-2000,10))
pl.legend(loc=(1.02,0.5))
pl.show()