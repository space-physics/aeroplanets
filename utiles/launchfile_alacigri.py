#!/usr/bin/env python
# -*- coding:Utf-8 -*-
from __future__ import division
from pylab import *
""" Main library for the pair production """
import os
import sys
import re
from threading import Thread
from datetime import datetime
montableau = []
monpilote = []
nombretotal = 0

launchxmlaero = os.path.join(os.path.dirname(__file__),"launchxmlaero.py")

def LaunchAero(filename, threadnumber):
	global montableau, monpilote
	print  "nombre de threads : ",threadnumber
	fi = open(filename,'r')
	for line in fi:
		if line.strip() != "":
			montableau.append(line)
			monpilote.append(0)
# L'execution consiste a chercher une valeur à 0 dans le tableau, a la faire passer à 1 et à executer l'equivalent de montableau.
# Ensuite, il suffit de mettre à 2 quand c'est fini
	for i in range(0,threadnumber):
		mythread=RunPlaneto(i)
		mythread.start()
	print "FINI?"
	print nombretotal

class RunPlaneto(Thread):
	def __init__(self,test):
		Thread.__init__(self)
		self.number=test

	def search(self):
		for i in range(0,len(monpilote)):
			if monpilote[i] == 0:
				return i
		return -1	
	def run(self): # Run est appelle par start pour les threads
		ok = True
		while ok:
			number=self.search()
			if number == -1:
				ok = False
			else:
				monpilote[number]=1
				print "le thread ",self.number,"  prend le nombre ",number
				launch = "nohup " + launchxmlaero +  "   " + montableau[number]
				print launch
				os.system(launch)
		#		os.system(os.path.join(os.path.dirname(__file__),"launchxmlaero.py") + "   " + montableau[number])
				global nombretotal
				nombretotal += 1

def ftostr(truc):
	return "%(str)f"%{"str":truc}

if __name__ == "__main__":
	print "salut"
	if len(sys.argv)!=3:
		print "Usage: launchfile_alacigri.py file threadnb"
		sys.exit()
	LaunchAero(sys.argv[1], int(sys.argv[2]))


