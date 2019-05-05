#!/usr/bin/env python
# -*- coding:Utf-8 -*-
from math import *
import os
import sys
import re
#from operator import itemgetter
""" Reads the resultat of the pair production simulations, and transform it into a coherent file """


nb = len(sys.argv)
if nb < 2:
    print "Vous devez donner le nom du fichier de resultat"
    sys.exit()
fileggname = sys.argv[1]
files = [i for i in os.listdir("./") if (os.path.isfile(i))]  # and re.match("launch_.*",i))]
print files
association = {}
for i in files:
    file = open(i, 'a+')
    count = 0
    for li in file:
        if (count != 0):
            nume = re.findall(r'\d+(?:\.\d+)?|NaN|Inf', li)
            association[float(nume[0])] = nume
        count += 1

meslignes = []
meslignes.append("# Energy per pair \n# Incident energy  \n# Mean production energy \n# elos1 \n#  elos2 \n# elos3 \n# ilos1 \n#  ilos2 \n# ilos3 \n# conservation \n#\n")

for i in sorted(association.keys()):
    truc = ""  # "%(ener)f\t"%{"ener":i}
    truc += '\t'.join(association[i])
    truc += "\n"
    meslignes.append(truc)
    print truc
    print i, "=>", '\t'.join(association[i])

print fileggname
file = open(fileggname, 'w')
file.writelines(meslignes)
file.close()
print "fin"
