#!/usr/bin/env python
# -*- coding:Utf-8 -*-
from __future__ import division
from pylab import *
import sys
from StringIO import StringIO
from glob import glob
from string import *
from scipy.stats import norm
from scipy.optimize import leastsq
try:
    import xml.etree.ElementTree as ET # in python >=2.5
except ImportError:
    try:
        import cElementTree as ET # effbot's C module
    except ImportError:
        try:
            import elementtree.ElementTree as ET # effbot's pure Python module
        except ImportError:
            try:
                import lxml.etree as ET # ElementTree API using libxml2
            except ImportError:
                import warnings
                warnings.warn("could not import ElementTree "
                              "(http://effbot.org/zone/element-index.htm)")
                # Or you might just want to raise an ImportError here.





def deriv_spl(x,y,d0,dn):
	""" Return the derivative at each point, for the spline
	x and y defines the function
	d0 derivative at the first point (>1E33 = 0)
	dn derivative for the last point
	"""
	if len(x)!=len(y):
		raise "fait Chier"
	dy=zeros(len(y))
	u=zeros(len(y))
	if(d0>1E33):
		dy[0]=0
		u[0]=0
	else:
		dy[0]=-0.5
		u[0]=(3./(x[1]-x[0]))*( (y[1]-y[0])/(x[1]-x[0]) -d0 )
	for i in range(1,len(y)-1):
		sig= (x[i]-x[i-1])/(x[i+1]-x[i-1])
		p= sig*dy[i-1]+2.
		dy[i]= (sig-1.)/p
		u[i]= (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1])
		u[i]=(6.*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p
	qn=0.
	if(dn>1E33):
	
		u[-1]=0
	else:
		qn=0.5
		u[-1]=(3./(x[-1]-x[-2]))*(dn- (y[-1]-y[-2])/(x[-1]-x[-2]) )

	dy[-1]=(u[-1]-qn*u[-2])/(qn*dy[-2]+1.)
	tmp=range(len(y)-1)
	tmp.reverse()
	print dy
	for i in tmp:
		dy[i]=dy[i]*dy[i+1]+u[i]
	print dy
	return dy



def splint(x,y,dy,newx):
	""" Interpolate the function x-y at newx. use the derivative computed by the precedent function """
	klo=0
	khi=len(y)-1

	while( khi-klo>1):
		k=(khi+klo) >>1
		if(x[k]<newx):
			klo=k
		else:
			khi=k
	h=x[khi]-x[klo]
	if(h==0):
		raise "et encore merde"
	a=(x[khi]-newx)/h
	b=(newx-x[klo])/h
	return a*y[klo]+b*y[khi]+((a**3-a)*dy[klo]+(b**3-b)*dy[khi])*h**2/6.




def spl_interp(x,y,d0,dn,newx):
	dy=deriv_spl(x,y,d0,dn)
	newy=zeros(len(newx))
	for i in range(len(newx)):
		newy[i]=splint(x,y,dy,newx[i])
	return newy
	#return dy


def spl_interpn(x,y,newx):
	d0=(y[1]-y[0])/(x[1]-x[0])
	dn=(y[-1]-y[-2])/(x[-1]-x[-2])
	dy=deriv_spl(x,y,d0,dn)
	newy=zeros(len(newx))
	for i in range(len(newx)):
		newy[i]=splint(x,y,dy,newx[i])
	return newy
	#return dy


def spl_interpexp(x,y,newx):
	return exp(spl_interpn(x[::-1],y[::-1],newx))


def GetArray(vNode,pname):
#	print vNode.find(pname).text
	return loadtxt(StringIO(vNode.find(pname).text.replace("\n"," ")))

def Getvalue(vNode,pname):
	return float(vNode.find(pname).text)



def ReadReport(repname):
	root = ET.parse(repname).getroot()
	chi2v = Getvalue(root,"chi2v")
	calibration = Getvalue(root,"calibration")
	print "Chi2 : ",chi2v,"Calibration",calibration
	alts = GetArray(root,"Specie/altitudes")
	vals = GetArray(root,"Specie/logvalues")
	return alts,vals, chi2v,calibration


def printhisto(data):
	hist(data,10)
	
def finderro(data):
	mu=data.mean()
	sigma=data.std()
#		mu,sigma=norm.fit(data)
	return mu,sigma

def plotnaltuncert(data,alt,nolog=False):
	newl = len(alt)

	mu = zeros((newl))
	sig = zeros((newl))

	for i in range(newl):
		mu[i],sig[i] = finderro(data[i,:])
	if not nolog:
		xscale("log")
	errorbar(mu,alt,xerr=sig)

def ExpProfile(pvals,z):
	ValMax=pvals[0]
	Texo=pvals[1]
	alt0=120
	R=3396.2
	mamu=44.001
	go=3.71
	SCALEH_CONST=8.3144727E5
	gamma=mamu*go*(1+alt0/R)**(-2)/(SCALEH_CONST*Texo)*1E5
	H=1/gamma
	return ValMax*exp(-(z-alt0)/H)
	
def ExpMin(pvals,z,ycompar):
	return (ycompar-ExpProfile(pvals,z))/ycompar

def expe(zalts,zmeasu):
	pretrieve=[1E11,300]
	pfinal=leastsq(ExpMin,pretrieve,args=(zalts,zmeasu))
#	print "Initial :",pretrieve
	print "Retrieved : ",pfinal
#	xscale("log")
#	plot(zmeasu,zalts,label="Data")
#	plot(ExpProfile(pfinal[0],zalts),zalts,label="Retrieved")
#	legend()
#	show()
	dens0,texo = pfinal[0]
	print dens0,texo
	return dens0,texo



def findrange(alt,alt0):
	""" returns the min, max range in altitude so that we are above alt0"""
	pos = -1
	if alt[-1]>alt[0]:
		rmax = len(alt)-1
		rmin = rmax
		while(rmin > -1 and alt[rmin]> alt0):
			rmin-=1
	else:
		rmin = 0
		rmax = rmin 
		while(rmax < len(alt) and alt[rmax] > alt0):
			rmax += 1
	return rmin,rmax



def ErrorTexo(data,alt,rmin,rmax):
	newl = len(data[0,:])
	dens0 = zeros((newl))
	Texo = zeros((newl))
	for i in range(newl):
		dens0[i],Texo[i] = expe(alt[rmin:rmax],data[rmin:rmax,i])
	print "Density uncertainty:", finderro(dens0)
	print "Texo uncertainty:", finderro(Texo)




if "__main__"==__name__:
	print "Salut les gars"
	files = glob("RapportFit.xml*")
	print files
	newalts = arange(80,200,1)
#	xscale("log")
	if not (len(files)>0):
		sys.exit()
	altitude, v, c2, ca = ReadReport(files[0])

	values = zeros((len(altitude),len(files)))
	newvalues = zeros((len(newalts),len(files)))
	chi = zeros(len(files))
	calibration = zeros(len(files))

	j = 0
	for i in files:
		print i
		alts,vals,chis, cals = ReadReport(i)
		values[:,j] = vals
		chi[j]=chis
		calibration[j]=cals
		newvalues[:,j]=spl_interpexp(alts,vals,newalts)
		#plot(spl_interpexp(alts,vals,newalts),newalts)
		j+=1
	
	rmin,rmax = findrange(newalts,135)
	print rmin,rmax,newalts[rmin],newalts[rmax]
	print newalts
	ErrorTexo(newvalues,newalts,rmin,rmax)
	#printhisto(calibration)
	#plotnaltuncert(newvalues,newalts)
	#plotnaltuncert(values,altitude,True)
	
	#print altitude
	#print values
	
	#show()
