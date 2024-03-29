#!/usr/bin/env python
# -*- coding:Utf-8 -*-
import os
import sys
import scipy
import re
from pylab import *
import pylab
from scipy.stats import norm


def string_to_tblnb(str):
    # Match a number  [-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?
    # return [float(i)  for i in re.findall(r'\d+(?:\.\d+)?',str)]
    #	print re.findall(r'([-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?)',str)
    return [float(i[0]) for i in re.findall(r'([-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?)', str)]


if len(sys.argv) == 1:
    print "zut", sys.argv
    print "you must define a regular expression to check your files"
    sys.exit(1)

regex_lect = sys.argv[1]


filenames = [i for i in os.listdir("./")
             if (os.path.isfile("./"+i))
             and re.match(regex_lect+"$", i)]

if len(filenames) == 0:
    print "no match for your file"
    sys.exit(1)


datas = []
altitude = array([])
defalt = False
skipped = 0
for i in filenames:
    print i

    fa = os.popen("grep -ri \"Electron impact ionization error\" "+i).read()
    value = string_to_tblnb(fa)[0]

    if(abs(value) > 15):
        skipped += 1
        print fa, " : SKIPPED"
        continue

    fi = os.popen("grep -ri nan "+i).read()
    if(len(fi) > 0):
        skipped += 1
        continue
    tmpdat = pylab.loadtxt("./"+i)
    if(not defalt):
        defalt = True
        altitude = tmpdat[:, 0]

    datas.append(tmpdat[:, 1:])

data = array(datas)

# print "bye"
fileintro = open("./"+filenames[0], 'r')
meslignes = []
for lines in fileintro:
    if(re.match("^#", lines)):
        meslignes.append(lines)


print data.shape
print altitude
(nbstat, nbalt, nbelem) = data.shape
sauvegarde = zeros((nbalt, nbelem*2+1))  # The number of elements corresponds to altitude+ elems+error
sauvegarde[:, 0] = altitude


for elem in range(nbelem):

    value = zeros(nbalt)
    error = zeros(nbalt)
    for alt in range(nbalt):
        tbla = data[:, alt, elem]
        mean = tbla.mean()
        std2 = tbla.std()*2
        #tbl=[tbla[i] for i in range(len(tbla)) if  abs(tbla[i]-mean)/mean<0.01   ]
        tbl = [tbla[i] for i in range(len(tbla)) if tbla[i] < mean+std2]
        print alt, len(tbla), len(tbl)
        mu, sigma = norm.fit(tbl)
        value[alt] = mu
        error[alt] = sigma
    print "position :", elem, 2*elem, 2*elem+1, "max :", nbelem*2
    sauvegarde[:, 2*elem+1] = value
    sauvegarde[:, 2*elem+2] = error

savetxt("montecarlo.dat", sauvegarde)
print "SKIPPED:", skipped

truc = file("montecarlo.dat", "a")
truc.writelines(meslignes)
truc.close()
