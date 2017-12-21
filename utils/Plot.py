#!/usr/bin/env python
from pylab import *
import pylab
from optparse import OptionParser
import sys
import re

parser = OptionParser()

xlog=False
ylog=False

lalegend=False

parser.add_option("-f","--file",dest="filename",help="The input file ",type="string")
parser.add_option("-o","--outfile",dest="imagename",help="The output image",type="string")
parser.add_option("--xlog",action="store_true",dest="xlog")
parser.add_option("--ylog",action="store_true",dest="ylog")
parser.add_option("--legend",action="store_true",dest="legend",help="If you want to see the legend")
parser.add_option("--invert",action="store_true",dest="invert",help="If you want to invert the x and y axis. Example : to plot the densities")
parser.add_option("--axis",dest="axis",help="Set the axis of your plot. Must be of the form \"[xmin,xmax,ymin,ymax]\"",type="string")





(options, args)=parser.parse_args()
#print options
#print "args"
#print args
if(options.legend):
	lalegend=True

filename=options.filename

if(not filename):
	if(len(args)):
		print("Dans les args : ")
		filename=args[0]
	else:
		print("Pas de fichiers")
		sys.exit()



file=open(filename,'r')

is_title=False
is_x=False
is_y=False
legende=[]



for lines in file:
	if(re.match("^#",lines)):
		if(not is_title):
			title(lines[2:])
			is_title=True
		else:
			if(not is_x):
				xlabel(lines[2:])
				is_x=True
			else:
				if(not is_y):
					ylabel(lines[2:])
					is_y=True
				else:
					legende.append(lines[2:])
	
file.close()



print("salut les gars",legende)




data=pylab.loadtxt(filename)



if options.invert:
	if not lalegend:
		for i in range(1,len(data[0])):
			plot(data[:,i],data[:,0])
	else:
		for i in range(1,len(data[0])):
			plot(data[:,i],data[:,0],label=legende[i-1])

else:
	if not lalegend:
		x=data[:,0]
		y=data[:,1:]
		plot(x,y)
	else:
		for i in range(1,len(data[0])):
			plot(data[:,0],data[:,i],label=legende[i-1])



if(options.xlog):
	xscale("log")

if(options.ylog):
	yscale("log")

if(options.axis):
	tableau=eval(options.axis)
	axis(tableau)

imagename=options.imagename
if(not imagename):
	imagename=filename+".png"
	if(lalegend):
		legend()
	savefig(imagename)
	print("show")
	show()
else:
	if(lalegend):
		legend()
	savefig(imagename)
	show()






