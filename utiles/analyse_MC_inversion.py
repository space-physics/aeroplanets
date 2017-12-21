#!/usr/bin/env python
# -*- coding:Utf-8 -*-
from glob import glob
from StringIO import StringIO
from pylab import *
from scipy.stats import norm
import os
import sys
import re
try:
	import xml.etree.ElementTree as ET # in python >=2.5
except ImportError:
	try:
		import cElementTree as ET 
	except ImportError:
		try:
			import elementtree.ElementTree as ET # effbot's pure Python module
		except ImportError:
			try:
				import lxml.etree as ET # ElementTree API using libxml2
			except ImportError:
				import warnings
				warnings.warn("could not import ElementTree "\
	                              "(http://effbot.org/zone/element-index.htm)")
                # Or you might just want to raise an ImportError here.


files = glob("report.xml*")
print files

mat = []
cal = []
chi = []
for i in files:
	print i
	handler = ET.parse(i)
	bal0=handler.find("/chi2v")
	chi.append(float(bal0.text))
	bal=handler.find("/Specie/logvalues")
	print bal.text
	data = loadtxt(StringIO(bal.text))
	bal2=handler.find("/calibration")
	calibration = float(bal2.text)
	mat.append(data)
	cal.append(calibration)

vals = matrix(mat)
cali = array(cal)
chiv = array(chi)
#print vals[:,0]

print "Points analysis"
for i in range(len(files)):
	print i,")------"
	ptx = vals[:,i]
	mu,sigma = norm.fit(ptx)
	print mu,sigma
	print ptx.mean(),ptx.std()

	print "--------------------"




print "Calibration analysis :"
mu,sigma = norm.fit(cali)
print cali.mean(),cali.std()
print mu,sigma

print "Chi2v analysis :"
mu,sigma = norm.fit(chiv)
print chiv.mean(),chiv.std()
print mu,sigma
print "Min", min(chiv),"Max",max(chiv)

print "bye"








