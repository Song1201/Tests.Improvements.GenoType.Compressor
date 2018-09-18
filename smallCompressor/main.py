#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 15:50:43 2018

@author: songchen
"""

import numpy as np
import vcf
import compressLib as cl
import importlib
import time

#%%
importlib.reload(cl)



#%%
cl.saveTemplate(filename='simulatedPopulation1.vcf.gz',template='template.vcf')

#%%
bitVectors,POSs,refs,alts,sampleNames = cl.readVcfgzToBitVectors(filename
                                        ='gene1000chr22.vcf.gz',startPos=16800001,
                                        endPos=17000000)
#%%
colBitVecs = cl.columnToRow(bitVectors)
bitVectors = colBitVecs
IDs = list(range(len(bitVectors))) # Initialize IDs
allHamsNext = np.array([bitVectors[ID].hamming_distance(bitVectors[ID+1]) for ID in
                  IDs[0:-1]],dtype=np.int16)
totalHam = np.sum(allHamsNext)
#%%
startTime = time.time()
heurisBitVecs,perm = cl.heurisHamming(bitVectors)
endTime = time.time()
print('Runtime:',endTime-startTime)
#%%
bitVectors, sampleDict = cl.countNonRef('gene1000chr22.vcf.gz',startPos=16800001,
                                        endPos=17000000)

#%%
colBitVecs = cl.columnToRow(bitVectors)

# How the number of Variant sites increase with respect to the number of samples
midBitVec = colBitVecs[0]
variantNo = np.array([])
sampleNo = np.array([])
for i in range(1,len(colBitVecs)):
    if i % 100 == 0:
        variantNo = np.append(variantNo,midBitVec.count_bits())
        sampleNo = np.append(sampleNo,i)
    midBitVec = midBitVec | colBitVecs[i]

variantNo = np.append(variantNo,midBitVec.count_bits())
sampleNo = np.append(sampleNo,len(colBitVecs))

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Define function for curve
def fExponential(n,a,b,c):
    return a*n**b + c 

def faLogN(n,a,b):
    return a*np.log(n) + b

popt, dump = curve_fit(fExponential,sampleNo,variantNo)
print(popt)
plt.plot(sampleNo,variantNo,'b-',label='Original Data')
plt.plot(sampleNo,fExponential(sampleNo,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of non-reference genotypes')
plt.ylabel('number of samples')
plt.savefig(fname='non_ref.eps',dpi=600,format='eps')
plt.show()
#%%
bitVectors = colBitVecs
IDs = list(range(amountBitVecs)) # Initialize IDs

startTime = time.time()
heurisBitVecs,perm = cl.heurisHamming(bitVectors)
endTime = time.time()
print('Runtime:',endTime-startTime)
#%%
allHamming = 0
for i in range(len(heurisBitVecs)-1):
    allHamming += heurisBitVecs[i].hamming_distance(heurisBitVecs[i+1])
print(allHamming)

#%%
allHamming = 0
for i in range(len(colBitVecs)-1):
    allHamming += colBitVecs[i].hamming_distance(colBitVecs[i+1])
print(allHamming)

#%%
unfitPos = np.load('unfitPos0.2_0.1.dat')
fitPos = np.load('fitPos0.2_0.1.dat')
totalHams = np.zeros(unfitPos.size+1)
for i in range(unfitPos.size):
    if i % 100 == 0:
        print('Running',i)
    totalHam = 0
    for j in range(len(colBitVecs)-1):
        totalHam += colBitVecs[j].hamming_distance(colBitVecs[j+1])
    totalHams[i] = totalHam
    colBitVecs[int(unfitPos[i])],colBitVecs[int(fitPos[i])] = colBitVecs[int(fitPos[i])].deep_copy(),colBitVecs[int(unfitPos[i])].deep_copy()
    
totalHam = 0
for j in range(len(colBitVecs)-1):
    totalHam += colBitVecs[j].hamming_distance(colBitVecs[j+1])
    
totalHams[unfitPos.size] = totalHam
totalHams.dump('totalHams0.2_0.1.dat')

#%%
IDs = list(range(5008))
for i in range(unfitPos.size):
    IDs[int(unfitPos[i])],IDs[int(fitPos[i])] = IDs[int(fitPos[i])], IDs[int(unfitPos[i])]

moved = 0    
for i in range(len(IDs)):
    if i != IDs[i]:
        moved += 1
print(moved)