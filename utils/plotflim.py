#!/usr/bin/env python
from pylab import *
import pylab
from optparse import OptionParser
import sys
import re
import os
import glob



def FileToPlot(filename,outname):
	clf()
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
	data=pylab.loadtxt(filename)
	for i in range(1,len(data[0])):
		plot(data[:,i],data[:,0],label=legende[i-1])
	xscale("log")
	legend()
	savefig(outname)


if len(sys.argv)==1:
	print "zut",sys.argv
	print "you must define a regular expression to check your files"
	sys.exit(1)


files=glob.glob(sys.argv[1])

print files

for i in files:
	number=re.findall(r"[0-9]+",i)[0]
	nb=number.zfill(5)
	FileToPlot(i,"tmp"+nb)
os.system("ffmpeg -r 2 -i tmp%5d.png test.ogv")
os.system("rm tmp[0-9]*.png")
#os.system("ffmpeg -r 10 -b 1800 -i tmp%05d test1800.mp4")


















