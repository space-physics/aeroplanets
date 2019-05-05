#!/usr/bin/env python
# -*- coding:Utf-8 -*-
from __future__ import division
from pylab import *
import os
import sys
import sqlite3 as lite


#############################
##
# To analyse the files
##
#############################


def returnAB(filename, wavea, waveb):
    """ Reads the filename, and finds the values of the different emissions in R, for the wavelength a and b
    return true and the two values if the system is successful
    return false and whatever has been found and/or 0 if not """
    data = loadtxt(filename)
    isaok = False
    isbok = False
    a = 0
    b = 0
    for i in range(len(data[:, 0])):
        if abs(wavea - data[i, 0]) < 1E-10:
            isaok = True
            a = data[i, 1]
        if abs(waveb - data[i, 0]) < 1E-10:
            isbok = True
            b = data[i, 1]

    return isaok & isbok, a, b


#############################
##
# create the database
##
#############################

def LfileOpenDB():
    """ Open the  a database"""
    if not os.path.isfile("superdb.db"):
        print "You should start the program with a negative option at the beginning to create the database"
        sys.exit()
        return False

    co = lite.connect("superdb.db")

    with co:
        cur = co.cursor()
        return co, cur
    return False


def LfileCreateDB(lfile):
    """ Creates a database"""
    co = lite.connect("superdb.db")
    with co:
        cur = co.cursor()
        dicolist = ReadArgFile(lfile)

        if len(dicolist) == 0:
            print "Your list is empty"
            sys.exit()
        keys = dicolist[0].keys()

        createstr = " FLOAT, ".join(keys)
        createstr += " FLOAT"

        # print createstr
        cur.execute("DROP TABLE IF EXISTS Simulations;")
        cur.execute("CREATE TABLE Simulations(Id INTEGER PRIMARY KEY, " + createstr + ")")
        # Orbit INT, Ls FLOAT, F107 FLOAT, Sza FLOAT, Alt FLOAT, Lat FLOAT, Long FLOAT, Localtime FLOAT, Band INT);")
        print "work on file"
        # Now, we work on the file
        for dico in dicolist:
            stk = "("
            stv = "("
            tmp = 0
            for ke in dico.keys():
                if tmp == 0:
                    stk += ke
                    stv += dico[ke]
                else:
                    stk += " , " + ke
                    stv += " , " + dico[ke]
                tmp += 1

            stk += ")"
            stv += ")"
            print stk
            print stv
            cur.execute("INSERT INTO Simulations" + stk + " VALUES " + stv + " ; ")
        co.commit()
        return co, cur
    return False


def ReadArgFile(lfile):
    dicolist = []
    f = open(lfile)
    for line in f:
        dico = {}
        args = line.split()
        dico["Key"] = args[1]
        for i in args[4:]:
            vals = i.split('=')
            dico[vals[0]] = vals[1]
        # print dico
        dicolist.append(dico)
    return dicolist

#######################################
###
# Analyse the database
###
#######################################


def ReadVariations(co, cur, key):
    """ Read the database and return a table with the variations in value of the parameter in key """
    cur.execute("SELECT " + key + " from Simulations;")
    co.commit()
    dic = {}
    while True:
        row = cur.fetchone()
        if row == None:
            break
        ro = row[0]
        dic[ro] = 0
    return dic.keys()


def ptype_study(co, cur, wa, wb):
    """ Plot the ratio a / b in function of the type"""
    # We study the variations of the main parameters
    ptypevar = ReadVariations(co, cur, "electron_ptype")
    print "electron_ptype", len(ptypevar), ptypevar
    eevar = ReadVariations(co, cur, "electron_energy")
    print "electron_energy", len(eevar), eevar
    E0var = ReadVariations(co, cur, "electron_E0")
    print "electron_E0", len(E0var), E0var
    n2var = ReadVariations(co, cur, "N2x")
    print "N2x", len(n2var), n2var

    for ty in ptypevar:
        print "SELECT (Key, electron_energy, electron_E0) FROM Simulations WHERE (N2x = 1 and electron_ptype = " + str(ty) + ") ORDER BY electron_E0;"
        cur.execute("SELECT Key , electron_energy , electron_E0  FROM Simulations WHERE (N2x = 1 and electron_ptype = " +
                    str(ty) + ")  ORDER BY electron_E0;")
        co.commit()
        E0 = []
        ratio = []
        while True:
            row = cur.fetchone()
            if row == None:
                break
            key = row[0]
            ezer = row[2]
            ret, a, b = returnAB("spectrum.out" + str(int(key)), wa, wb)
            if ret:
                E0.append(ezer)
                ratio.append(a / b)
    #	dic[ty] = (E0, ratio)
        plot(E0, ratio, label="Precip type :" + str(ty))
    legend()
    show()


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print "Usage: launch_analysis.py launchfilename wavelength_a wavelength_b option"
        print "You should launch that file in the outputdir"
        print "The launchfilename is the name of the file used for creating the simulations"
        print "The wavelength _a and _b are the wavelength with which we will make a ratio _a / _b"
        print "The option gives the kind of analysis made (against energy, against E0...)"
        # Nb: we will have the possibility to use the MC factor to make a different analysis...
        # The fact that this factor exists also ensure some funny part: we will be able to
        # perform the division for some ratio of the same computation
        sys.exit()
    lfile = sys.argv[1]
    wa = float(sys.argv[2])
    wb = float(sys.argv[3])
    option = int(sys.argv[4])

    # ReadArgFile(lfile)
    # sys.exit()
    if option <= 0:
        co, cur = LfileCreateDB(lfile)
    else:
        co, cur = LfileOpenDB()
    option = abs(option)
    # We can open the database directly

    ptype_study(co, cur, wa, wb)
