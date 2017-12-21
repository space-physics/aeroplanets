#!/usr/bin/env python
from matplotlib.pyplot import figure, show
import re
import numpy as np
from argparse import ArgumentParser

xlog=False
ylog=False

lalegend=False

p = ArgumentParser()
p.add_argument("infn",help="The input file.")
p.add_argument("-o","--outfn",help="The output image")
p.add_argument("--xlog",action="store_true")
p.add_argument("--ylog",action="store_true")
p.add_argument("--legend",action="store_true",help="If you want to see the legend")
p.add_argument("--invert",action="store_true",help="If you want to invert the x and y axis. Example : to plot the densities")
p.add_argument("--axis",help="Set the axis of your plot. Must be of the form [xmin,xmax,ymin,ymax]")
p = p.parse_args()

# %%
is_title=False
is_x=False
is_y=False
legende=[]

fig = figure()
ax = fig.gca()


with open(p.infn,'r') as f:

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
data=np.loadtxt(p.infn)

if p.invert:
	if not p.legend:
		for i in range(1,len(data[0])):
			ax.plot(data[:,i],data[:,0])
	else:
		for i in range(1,len(data[0])):
			ax.plot(data[:,i],data[:,0],label=legende[i-1])

else:
	if not p.legend:
		x=data[:,0]
		y=data[:,1:]
		ax.plot(x,y)
	else:
		for i in range(1,len(data[0])):
			ax.plot(data[:,0],data[:,i],label=legende[i-1])



if p.xlog:
	ax.set_xscale("log")

if p.ylog:
	ax.set_yscale("log")

if p.axis:
	tableau=eval(p.axis)
	ax.axis(tableau)

if p.legend:
   ax.legend(loc='best')

# %%
if not p.outfn:
	imagename=p.infn+".png"

fig.savefig(imagename)


show()







