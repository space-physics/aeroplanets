#!/usr/bin/env python
import scipy.stats
import matplotlib.pyplot as plt
from optparse import OptionParser
#from powerlaw import powerlaw
import xml.etree.ElementTree as ET # in python >=2.5
import numpy as np
from io import StringIO
import sys


def PlotShiraiNode(vNode):
	leg=""
	if 0==len(vNode.findall("legend")):
		leg=vNode.attrib["name"]
	else:
		leg=vNode.find("legend").text
	try:
		threshold=float(vNode.attrib["threshold"])
	except:
		threshold=0
	uncertainty=float(vNode.find("uncertainty").text)
	Emin=float(vNode.find("Emin").text)
	Emax=float(vNode.find("Emax").text)

	tid=vNode.find("Equation").attrib["article_id"]
	tip=int(vNode.find("Equation").attrib["type"])
	params=np.loadtxt(StringIO(vNode.find("params").text.replace("\n"," ")))

#	def __init__(self,emin,emax,threshold,eqtype,dataeq):
	if(tid=="CH4"):
		shirai=NewShiraiCH4(Emin,Emax,threshold,tip,params)
	if(tid=="CO2"):
		shirai=NewShiraiCO2(Emin,Emax,threshold,tip,params)
	if(tid=="N2"):
		shirai=NewShiraiN2(Emin,Emax,threshold,tip,params)

	ene= np.arange(Emin,Emax,(Emax-Emin)/100.)
	ene= np.array(scipy.stats.powerlaw(Emin,Emax,100))
	print(type(ene))
	cross=shirai.ReturnCrs(ene*1E-3)
	datauncert=cross*uncertainty/100.
	plt.errorbar(ene,cross,yerr=datauncert,label=leg)


def PlotStdNode(vNode):
	leg=""
	if(0==len(vNode.findall("legend"))):
		try:
			leg=vNode.attrib["name"]
		except:
			leg="Elastic"
	else:
		leg=vNode.find("legend").text

	fact=1
	if("fact" in list(vNode.find("Egrid").keys())):
		fact=float(vNode.find("Egrid").attrib.get("fact"))
#	print "Votre facteur :",fact
#	print loadtxt(StringIO((vNode.find("Egrid").text).replace("\n"," ")))
	dataenergy=np.loadtxt(StringIO(vNode.find("Egrid").text.replace("\n"," ")))*fact
	fact=1
	if("fact" in list(vNode.find("Cross").keys())):
		fact=float(vNode.find("Cross").attrib.get("fact"))
	print("Votre facteur :",fact)
	datacrs=np.loadtxt(StringIO(vNode.find("Cross").text.replace("\n"," ")))*fact
#	print datacrs

	uncertainty=0
	datauncert=np.zeros((len(datacrs)))
	if("uncertainty" in list(vNode.find("Cross").keys())):
		uncertainty=vNode.find("Cross").attrib.get("uncertainty")
		if (uncertainty.find("%")):
			uncert=float(uncertainty.replace("%",""))/100.
			print("Facteur d'incertitude",uncert)
			datauncert=datacrs*uncert
		else:
			value=float(uncertainty)*fact
			print("Valeur d'incertitude")
			datauncert=np.ones((len(datacrs)))*value

#	print datauncert

#	errorbar(datacrs,dataenergy,yerr=datauncert,label=leg)
	plt.errorbar(dataenergy,datacrs,yerr=datauncert,label=leg)



def CheckNode(vNode,nodename=""):
    if nodename=="":
        nodename=vNode.attrib["name"]
    if("uncertainty" in list(vNode.find("Cross").keys())):
        return True
    if("fact_uncertainty" in list(vNode.find("Cross").keys())):
        return True

    raise RuntimeError("Problem for the node ",nodename)


def CheckShirai(vNode,name=""):
    if(len(vNode.findall("uncertainty"))!=1):
        if name=="":
            e = "Problem with {}".format(vNode.attrib["name"])
        else:
            e = "Problem with {}".format(name)

        raise RuntimeError(e)


def ActiveNode(vNode):
	vNode.find("Cross").attrib["setMC"]="active"

def DeActiveNode(vNode):
	vNode.find("Cross").attrib["setMC"]="inactive"

def ActiveShirai(vNode):
	if(not vNode.find("SetMCActive")):
		ET.SubElement(vNode,"SetMCActive")

def DeActiveShirai(vNode):
	for i in (vNode.findall("SetMCActive")):
		vNode.remove(i)


def ForceNode(vNode,st):
	if("uncertainty" in list(vNode.find("Cross").keys())):
		return
	if("fact_uncertainty" in list(vNode.find("Cross").keys())):
		return
	vNode.find("Cross").attrib["uncertainty"]=st
def ForceShirai(vNode,st):
	if(len(vNode.findall("uncertainty"))!=1):
		unc=ET.SubElement(vNode,"uncertainty")
		unc.text=st

def ToFactor(vNode):
	uncert=vNode.find("Cross").attrib["uncertainty"]
	del vNode.find("Cross").attrib["uncertainty"]
	vNode.find("Cross").attrib["fact_uncertainty"]=uncert

def ToUncert(vNode):
	uncert=vNode.find("Cross").attrib["fact_uncertainty"]
	del vNode.find("Cross").attrib["fact_uncertainty"]
	vNode.find("Cross").attrib["uncertainty"]=uncert



if __name__=="__main__":
	if(len(sys.argv)<2):
		print("veuillez donner un nom de fichier")
		sys.exit()


	parser = OptionParser()

	parser.add_option("-f","--file",dest="filename",help="The input file ",type="string")
	parser.add_option("-c","--check-mc",help="Check the existence of MonteCarlo",action="store_true",dest="check")
	parser.add_option("-a","--activate",help="Activate the Monte Carlo system",action="store_true",dest="active")
	parser.add_option("-d","--deactivate",help="DeActivate the Monte Carlo system",action="store_true",dest="deactive")
	parser.add_option("--force",help="Force an uncertainty in %",dest="force",type="string")
	parser.add_option("--tofactor",help="The uncertainties are converted to factor uncertainties",action="store_true",dest="tofactor")
	parser.add_option("--touncert",help="The factor uncertainties are converted to standard uncertainties",action="store_true",dest="touncert")
	parser.add_option("-o","--outfile",dest="outname",help="The output image",type="string")

	(options, args)=parser.parse_args()

	filename=options.filename

	if(not filename):
		if(len(args)):
			print("Dans les args : a")
			filename=args[0]
		else:
			print("Pas de fichiers")
			sys.exit()
	print("Fichier",filename)


	check=options.check
	active=options.active
	deactive=options.deactive

	touncert=options.touncert
	tofactor=options.tofactor
	if tofactor and touncert:
		raise ValueError("you cannot convert to factor and to uncert at the same time")
	if active:
		check=True
		if(not options.outname):
			raise ValueError("Error: out file not set")
		if(deactive):
			raise ValueError("You cannot activate and deactivate at the same time (dumb)")
	force=False
	forcesh=""
	forcen=""
	if options.force:
		force=True
		forcesh=(options.force).strip()
		forcen=forcesh+"%"
		print("Force",forcesh,forcen)



	tree=ET.parse(filename)
	root=tree.getroot()
	processlist=root.findall(".//Process")
	print("Nous avons trouve ",len(processlist),"processus")

	for proc in processlist:
		if(0==len(proc.findall("Shirai"))):
			if tofactor:
				ToFactor(proc)
			if touncert:
				ToUncert(proc)
			if force:
				ForceNode(proc,forcen)
			if check:
				CheckNode(proc)
			if active:
				ActiveNode(proc)
			if deactive:
				DeActiveNode(proc)
		else:
			if force:
				ForceShirai(proc,forcesh)
			if check:
				CheckShirai(proc)
			if active:
				ActiveShirai(proc)
			if deactive:
				DeActiveShirai(proc)

	processlist2=root.findall(".//ElasticCrs")
	for proc in processlist2:
		if(0==len(proc.findall("Shirai"))):
			if tofactor:
				ToFactor(proc)
			if touncert:
				ToUncert(proc)
			if force:
				ForceNode(proc,forcen)
			if check:
				CheckNode(proc,"ElasticCrs")
			if active:
				ActiveNode(proc)
			if deactive:
				DeActiveNode(proc)
		else:
			if force:
				ForceShirai(proc,forcesh)
			if check:
				CheckShirai(proc,"ElasticCrs")
			if active:
				ActiveShirai(proc)
			if deactive:
				DeActiveShirai(proc)

	processlist3=root.findall(".//TotalCrs")
	for proc in processlist3:
		if(0==len(proc.findall("Shirai"))):
			if tofactor:
				ToFactor(proc)
			if touncert:
				ToUncert(proc)
			if force:
				ForceNode(proc,forcen)
			if check:
				CheckNode(proc,"TotalCrs")
			if active:
				ActiveNode(proc)
			if deactive:
				DeActiveNode(proc)
		else:
			if force:
				ForceShirai(proc,forcesh)
			if check:
				CheckShirai(proc,"TotalCrs")
			if active:
				ActiveShirai(proc)
			if deactive:
				DeActiveShirai(proc)


	if(options.outname):
		tree.write(options.outname)


