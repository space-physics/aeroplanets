#!/usr/bin/env python
# -*- coding:Utf-8 -*-
from __future__ import division
import re
import os
import sys
try:
    import xml.etree.ElementTree as ET  # in python >=2.5
except ImportError:
    try:
        import cElementTree as ET  # effbot's C module
    except ImportError:
        try:
            import elementtree.ElementTree as ET  # effbot's pure Python module
        except ImportError:
            try:
                import lxml.etree as ET  # ElementTree API using libxml2
            except ImportError:
                import warnings
                warnings.warn("could not import ElementTree "
                              "(http://effbot.org/zone/element-index.htm)")
                # Or you might just want to raise an ImportError here.


if len(sys.argv) < 4:
    sys.exit()

# 1 We read the parameters
args = sys.argv[4:]
filelog = sys.argv[1]  # Not used here, but necessary for CIGRI
template = sys.argv[3]  # The template to modify
nb = sys.argv[2]  # The number of the file (to simplify)
dico = {}
for i in args:
    vals = i.split('=')
    dico[vals[0]] = vals[1]

# We open the file
print "Template", template
handler = ET.parse(template)
root = handler  # .find("aero_main")
name = template + "%08i" % int(nb)
#tree = ET.ElementTree(root)
# tree.write(name)
# root.write(name)
# print "merde1", name
# sys.exit()
# print "merde"


extra = "# Init\n # " + " ".join(sys.argv) + "#Used\n# "
# 3 We look at the different modifications to perform in the file

# For IRI and MSIS and SUN : f107
if "f107" in dico.keys():
    b1 = handler.find("./sun/model/f107")
    b1.text = dico["f107"]
    extra += "f107_" + dico["f107"]

if "Ap" in dico.keys():
    b1 = handler.find("./planet/Ap")
    b1.text = dico["Ap"]
    extra += "Ap_" + dico["Ap"]


# Atmosphere system
if "N2x" in dico.keys():
    b1 = handler.find("./atmosphere/neutral/multiplicator/N2")
    b1.text = dico["N2x"]
    extra += "N2x_" + dico["N2x"]

if "O2x" in dico.keys():
    b1 = handler.find("./atmosphere/neutral/multiplicator/O2")
    b1.text = dico["O2x"]
    extra += "O2x_" + dico["O2x"]

if "Ox" in dico.keys():
    b1 = handler.find("./atmosphere/neutral/multiplicator/O")
    b1.text = dico["Ox"]
    extra += "Ox_" + dico["Ox"]

for sp in ["O", "O2", "N2"]:
    spdir = "./atmosphere/neutral/add_parametrized_species/" + sp + "/Model/"

    bwalt = "BW_" + sp + "_alt"
    if bwalt in dico.keys():
        print spdir
        b1 = handler.find(spdir + "alt")
        b1.text = dico[bwalt]
        extra += bwalt + "_" + dico[bwalt]

    bw = "BW_" + sp + "_dens"
    if bw in dico.keys():
        b1 = handler.find(spdir + "dens")
        b1.text = dico[bw]
        extra += bw + "_" + dico[bw]

    bw = "BW_" + sp + "_Texo"
    if bw in dico.keys():
        b1 = handler.find(spdir + "Texo")
        b1.text = dico[bw]
        extra += bw + "_" + dico[bw]

    bw = "BW_" + sp + "_To"
    if bw in dico.keys():
        b1 = handler.find(spdir + "To")
        b1.text = dico[bw]
        extra += bw + "_" + dico[bw]

    bw = "BW_" + sp + "_Shape"
    if bw in dico.keys():
        b1 = handler.find(spdir + "Shape")
        b1.text = dico[bw]
        extra += bw + "_" + dico[bw]


# Electron system

if "electron_used" in dico.keys():
    b1 = handler.find("./electron")
    if dico["electron_used"] != "0":
        print "We use the electrons"
        if not ET.iselement(handler.find("./electron/use_electron")):
            b1.append(ET.Element("use_electron"))
    else:
        print "suppress"
        try:
            for b2 in handler.findall("./electron/use_electron"):
                print b2
                b1.remove(b2)
        except:
            print "Nothing to remove"

if "electron_ptype" in dico.keys():
    b1 = handler.find("./electron/precipitation/use_model")
    b1.attrib["type"] = dico["electron_ptype"]
    extra += "electronptype_" + dico["electron_ptype"]

if "electron_energy" in dico.keys():
    b1 = handler.find("./electron/precipitation/use_model/entot")
    b1.text = dico["electron_energy"]
    extra += "electronenergy_" + dico["electron_energy"]

if "electron_E0" in dico.keys():
    b1 = handler.find("./electron/precipitation/use_model/E0")
    b1.text = dico["electron_E0"]
    extra += "electronE0_" + dico["electron_E0"]

if "electron_powlaw" in dico.keys():
    b1 = handler.find("./electron/precipitation/use_model/powlaw")
    b1.text = dico["electron_powlaw"]
    extra += "electronpowlaw_" + dico["electron_powlaw"]

if "electron_isotro" in dico.keys():
    b1 = handler.find("./electron/precipitation/use_model/isotro")
    b1.text = dico["electron_isotro"]
    extra += "electronisotro_" + dico["electron_isotro"]


# Proton system
if "proton_used" in dico.keys():
    b1 = handler.find("./proton")
    if dico["proton_used"] != "0":
        print "We use the protons"
        if not ET.iselement(handler.find("./proton/use_proton")):
            b1.append(ET.Element("use_proton"))
    else:
        print "suppress"
        try:
            for b2 in handler.findall("./proton/use_proton"):
                print b2
                b1.remove(b2)
        except:
            print "Nothing to remove"
    #b1 = handler.find("/proton/model/continuum/Egrid")

if "proton_ptype" in dico.keys():
    b1 = handler.find("./proton/proton_precip/precipitation/use_model")
    b1.attrib["type"] = dico["proton_ptype"]
    extra += "protonptype_" + dico["proton_ptype"]

if "proton_energy" in dico.keys():
    b1 = handler.find("./proton/proton_precip/precipitation/use_model/entot")
    b1.text = dico["proton_energy"]
    extra += "protonenergy_" + dico["proton_energy"]

if "proton_E0" in dico.keys():
    b1 = handler.find("./proton/proton_precip/precipitation/use_model/E0")
    b1.text = dico["proton_E0"]
    extra += "protonE0_" + dico["proton_E0"]

if "proton_powlaw" in dico.keys():
    b1 = handler.find("./proton/proton_precip/precipitation/use_model/powlaw")
    b1.text = dico["proton_powlaw"]
    extra += "protonpowlaw_" + dico["proton_powlaw"]

if "proton_isotro" in dico.keys():
    b1 = handler.find("./proton/proton_precip/precipitation/use_model/isotro")
    b1.text = dico["proton_isotro"]
    extra += "protonisotro_" + dico["proton_isotro"]


# 4 We write the extra info

extra += "#\n# "


b4 = handler.find("./ExtraInfo")
try:
    b4.text = extra
except:
    b3 = handler.find("./")
    b3.append(ET.Element("ExtraInfo"))
    b3.text = extra


# 5 We correct the outputs


def txtoutcorr(dat, outdir):
    (odir, fname) = os.path.split(dat)
    return os.path.join(outdir, fname)


def mkoutdir(dirname):
    try:
        os.makedirs(dirname)
    except:
        print "Dir already created"


if "outdir" in dico.keys():
    outdir = dico["outdir"]
    mkoutdir(outdir)
    corrbal = ["./chem/chem_atmo",
               "./emissions/emit_list",
               "./emissions/spectrum",
               "./emissions/emission/column_file",
               "./emissions/emission/profile_file",
               "./emissions/emission/extra_file",
               "./output/neutral_atmo",
               "./output/energy_conservation"
               ]

    for bal in corrbal:
        for b1 in handler.findall(bal):
            b1.text = txtoutcorr(b1.text, outdir)
# 6 We write the file

#tree = ET.ElementTree(root)
# tree.write(name)
root.write(name)
# 7 We launch the system
os.system("aero1d " + name + " " + nb)
# exemple
# launchxmlaero.py DayEarthProton.xml 0 proton_used=1 proton_energy=42
# launchxmlaero.py DayEarthProton.xml 1 proton_used=0 proton_energy=43
# launchxmlaero.py DayEarthProton.xml 2 proton_energy=44
