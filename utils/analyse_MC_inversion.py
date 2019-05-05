#!/usr/bin/env python
# -*- coding:Utf-8 -*-
from glob import glob
from io import StringIO
import numpy as np
from scipy.stats import norm
import xml.etree.ElementTree as ET

pat = "report.xml*"

files = glob(pat)
if not files:
    raise FileNotFoundError('did not find files with glob {}'.format(pat))

mat = []
cal = []
chi = []
for i in files:
    print(i)
    handler = ET.parse(i)
    bal0 = handler.find("/chi2v")
    chi.append(float(bal0.text))
    bal = handler.find("/Specie/logvalues")
    print(bal.text)
    data = np.loadtxt(StringIO(bal.text))
    bal2 = handler.find("/calibration")
    calibration = float(bal2.text)
    mat.append(data)
    cal.append(calibration)

vals = np.matrix(mat)
cali = np.array(cal)
chiv = np.array(chi)
# print vals[:,0]

print("Points analysis")
for i in range(len(files)):
    print(i, ")------")
    ptx = vals[:, i]
    mu, sigma = norm.fit(ptx)
    print(mu, sigma)
    print(ptx.mean(), ptx.std())

    print("--------------------")


print("Calibration analysis :")
mu, sigma = norm.fit(cali)
print(cali.mean(), cali.std())
print(mu, sigma)

print("Chi2v analysis :")
mu, sigma = norm.fit(chiv)
print(chiv.mean(), chiv.std())
print(mu, sigma)
print("Min", min(chiv), "Max", max(chiv))
