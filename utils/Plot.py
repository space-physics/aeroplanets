#!/usr/bin/env python
from matplotlib.pyplot import figure, show
import re
import numpy as np
from argparse import ArgumentParser
#
xlog=False
ylog=False
lalegend=False
#
p = ArgumentParser()
p.add_argument("path", help="The aero1d output directory to plot from")
p.add_argument("-o","--outfn", help="The output image")
p = p.parse_args()
# %%
is_title=False
is_x=False
is_y=False
legende=[]

fig = figure()
ax = fig.gca()

with open(p.path, 'r') as f:

    for lines in f:
	    if(re.match("^#",lines)):
		    if not is_title:
			    ax.set_title(lines[2:])
			    is_title=True
		    else:
			    if not is_x:
				    ax.set_xlabel(lines[2:])
				    is_x=True
			    else:
				    if not is_y:
					    ax.set_ylabel(lines[2:])
					    is_y=True
				    else:
					    legende.append(lines[2:])

print(legende)
# %%
data=np.loadtxt(p.path)

for i in range(1, len(data[0])):
	ax.plot(data[:,i], data[:,0], label=legende[i-1])

ax.set_xscale("log")
ax.set_xlim(1e-2,None)
ax.legend(loc='best')
# %%
if not p.outfn:
	imagename = p.path+".png"

fig.savefig(imagename)

show()
