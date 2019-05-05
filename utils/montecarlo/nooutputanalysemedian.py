#!/usr/bin/env python
# -*- coding:Utf-8 -*-
import os
import sys
import scipy
import re
from pylab import *
import pylab


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

    if(abs(value) > 5):
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
sauvegarde = zeros((nbalt, nbelem*5+1))  # The number of elements corresponds to altitude+ elems+error
sauvegarde[:, 0] = altitude


for elem in range(nbelem):

    value = zeros(nbalt)
    error1 = zeros(nbalt)
    error2 = zeros(nbalt)
    vmin = zeros(nbalt)
    vmax = zeros(nbalt)

    for alt in range(nbalt):
        tbl = data[:, alt, elem]
        value[alt] = median(tbl)
        # error[alt]=tbl.std()*2
        mask1 = (tbl < value[alt])
        error1[alt] = tbl[mask1].std()*2
        mask2 = (tbl > value[alt])
        error2[alt] = tbl[mask2].std()*2
        vmin[alt] = tbl.min()
        vmax[alt] = tbl.max()

    print "position :", elem, 5*elem, 5*elem+1, "max :", nbelem*2
    sauvegarde[:, 5*elem+1] = value
    sauvegarde[:, 5*elem+2] = error1
    sauvegarde[:, 5*elem+3] = error2
    sauvegarde[:, 5*elem+4] = vmin
    sauvegarde[:, 5*elem+5] = vmax

savetxt("montecarlomedian.dat", sauvegarde)
print "SKIPPED:", skipped

truc = file("montecarlomedian.dat", "a")
truc.writelines(meslignes)
truc.close()


xscale("log")
yscale("log")
print altitude[156]
tbl = data[:, 156, 0]
print "mean :", tbl.mean()
med = median(tbl)
print "median : ", med
# error[alt]=tbl.std()*2
mask1 = (tbl < med)
print "Error min ", tbl[mask1].std()*2
mask2 = (tbl > med)
print "Error max ", tbl[mask2].std()*2
print "min : ", tbl.min()
print "max : ", tbl.max()
hist(tbl, 1000)
show()
