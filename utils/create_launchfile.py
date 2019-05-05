#!/usr/bin/env python
from powerlaw import powerlaw
import sys
dico = {}

# The IRI, MSIS parameters
dico["f107"] = [100]
dico["Ap"] = [10]

# The atmosphere multiplication factors
dico["N2x"] = [1, 2, 3]
dico["O2x"] = [1]
dico["Ox"] = [1]

# The BatesWalker profiles parameters
#

# The electron system

dico["electron_used"] = [1]
dico["electron_ptype"] = [1]
dico["electron_energy"] = powerlaw(10, 10000, 50)  # arange(1000, 10000,250)
dico["electron_E0"] = powerlaw(10, 10000, 50)  # arange(1000, 10000, 250)
dico["electron_powlaw"] = [0]
dico["electron_isotro"] = [1]

# The Proton system

dico["proton_used"] = [0]
dico["proton_ptype"] = [1]
dico["proton_energy"] = [0]  # arange(1000, 10000,250)
dico["proton_E0"] = [0]  # arange(1000, 10000, 250)
dico["proton_powlaw"] = [0]
dico["proton_isotro"] = [1]


# All of that is multiplied by the number of monte carlo runs
dico["MCRun"] = [0]


# We create the table with all the parameters

def RecCreateFiles(dic, keys, pos):
    """ Returns an array with all the files"""
    if (pos == len(keys)):
        return []
    downarray = RecCreateFiles(dic, keys, pos + 1)
    myarr = []
    st = keys[pos]
    for val in dic[st]:
        print(st, val)
        if downarray == []:
            myarr.append(st + "=" + str(val))
        else:
            for dat in downarray:
                myarr.append(st + "=" + str(val) + "   " + dat)

    return myarr


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: create_launchfile.py templatename outdirname launchfilename")
        raise SystemExit()
    template = sys.argv[1]
    outdir = sys.argv[2]
    launchfile = sys.argv[3]
    myvallist = RecCreateFiles(dico, list(dico.keys()), 0)

    with open(launchfile, "w") as f:
        for i in range(len(myvallist)):
            print("Log.txt"+str(i), i, template.strip(), "outdir=" + outdir.strip(), myvallist[i], file=f)


# print dico.keys()
#outdir= "  superout   "
# print myvallist
# for i in range(len(myvallist)):
#	print "Log.txt"+str(i), i, "testtemplate.xml", "outdir=" + outdir.strip(), myvallist[i]
