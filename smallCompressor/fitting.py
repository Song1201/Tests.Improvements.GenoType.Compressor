# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 23:35:39 2018

@author: lison
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#%% Import data
decompressTime = np.array([1.059,2.484,4.476,10.56,33.98	,72.58,122.7,194.5,273.3,370,481.7,
                           609.7,772.5])
query23239Time = np.array([0.5268,0.5368,0.653,0.6282,0.9504,1.463,2.172,3.185,4.307,5.61,
                           7.104,8.926,11.35,])
query1mTime = np.array([0.7311,1.412,2.243,5.081,15.96,33.36,57.71,89.07,133.6,174.4,                  
                        227.6,287.6,363.4])
query1variantTime = np.array([0.5275,0.5248,0.5429,0.5657,0.5964,0.6111,0.6819,0.7639,0.8124,	
                              0.8281,0.8911,0.9244,1.037])
query1sampleTime = np.array([0.5539,0.5615,0.637,0.7946,1.087,1.444,1.894,3.041,3.374,
                             4.647,6.7,7.055,7.742])
sampleNumber = np.array([1,2,3,5,10,15,20,25,30,35,40,45,50])*10000
variantNumber = np.array([46825,62285,74396,94984,143562,189727,234773,278649,321775,363445,
                          403914,443394,481466])
record1m = np.array([18844,25670,31273,40851,64016,86097,107554,128427,148863,168563,187821,
                     206522,224651])
phase1 = np.array([1.05,0.59,0.66,0.61,0.59,0.54,0.7,0.44,0.56,0.62,0.6563,0.5,0.5])
phase2 = np.array([16.4,32.97,30.17,24.7,48.91,79.04,98.8,122.84,198.64,231.46,310.8,390.5,            
                   390.1])
phase3 = np.array([1.37,13.61,28.71,78.7,323.25,762.31,1380.1,2133.46,3039.6,4179.66,5347,
                   6878,8666])
phase4 = np.array([0.51,1.92,3.1,6.93,25.48,52.17,96.04,158.48,219.3,329.18,458.6,614.9,768])
phase5 = np.array([0.04,0.03,0.16,0.39,0.72,2.98,2.12,3.72,9.66,12.82,6.375,9.594,12.97])
permutation23239 = np.array([0.5625,9.25,24.09,69.56,279,643.8,1179,1812,2583,3601,4630,5981,7540])
real = np.array([426.7,484.3,522.6,609.8,935.9,1443,2165,3043,4087,5502,6893,8773,10911])
iteration = np.array([150,267,392,647,1292,1954,2661,3269,3891,4463,5066,5676,6248])
permPerIter = permutation23239/iteration
bitVectorNumber = np.array([93958,125356,150176,192726,294608,393194,491166])
perm1m = np.array([0.07281,0.4684,1.081,2.913,11.32,24.69,43.77,69.29,105.1,134.5, 
                          173.9,220.3,275.1])
perm1mPerRecord = perm1m/record1m
real1m = np.array([0.6732,1.292,2.146,4.713,15.83,33.02,57.41,89.79,136.5,176.5,226.6,289.6,
                   362.2])
real1mPerRecord = real1m/record1m
perm1mNew = np.array([0.073,0.4042,0.9201,2.46,9.971,21.21,37.62,57.68,87.83,116.4,152.9,196.1,
                      238.7])
real1mNew = np.array([0.6962,1.234,2.057,4.273,14.54,29.72,51.45,78.63,119.1,158.2,207.8,265,
                      324.4])
perm1mPerRecord = perm1m/record1m
perm1mNewPerRecord = perm1mNew/record1m
compresSize = np.array([2,4.1,6.2,12,31,58.5,89,128,178,232,293,358,448.7])

permutation = np.array([0.5625,9.25,24.09,69.56,279,643.8,1179,1812,2583,3601,4630,5981,7540])
real = np.array([426.7,484.3,522.6,609.8,935.9,1443,2165,3043,4087,5502,6893,8773,10911])
iteration = np.array([150,267,392,647,1292,1954,2661,3269,3891,4463,5066,5676,6248])
permPerIter = permutation/iteration
amountData = sampleNumber*variantNumber/1000000
queryPerRecord = decompressTime/variantNumber
real1mNewPerRecord = real1mNew/record1m

rest1mPerRecord = real1mPerRecord - perm1mPerRecord

query1mPerRecord = query1mTime/record1m


#%% Curve Fitting
# Define function for curve
def fExponential(n,a,b,c):
    return a*n**b + c 

def faLogN(n,a,b):
    return a*np.log(n) + b


#%%
popt, dump = curve_fit(fExponential,sampleNumber,decompressTime)
print(popt)
plt.plot(sampleNumber,decompressTime,'b-',label='Original Data')
plt.plot(sampleNumber,fExponential(sampleNumber,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('time(s)')
plt.savefig(fname='decompress.eps',dpi=600,format='eps')
plt.show()

#%%
popt, dump = curve_fit(fExponential,sampleNumber,query23239Time)
print(popt)
plt.plot(sampleNumber,query23239Time,'b-',label='Original Data')
plt.plot(sampleNumber,fExponential(sampleNumber,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('time(s)')
plt.savefig(fname='query23239.eps',dpi=600,format='eps')
plt.show()

#%%
popt, dump = curve_fit(fExponential,amountData,query1mTime)
print(popt)
plt.plot(amountData,query1mTime,'b-',label='Original Data')
plt.plot(amountData,fExponential(amountData,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('amount of data(MB)')
plt.ylabel('time(s)')
plt.savefig(fname='query1m.eps',dpi=600,format='eps')
plt.show()

#%%
popt, dump = curve_fit(fExponential,sampleNumber,amountData)
print(popt)
plt.plot(sampleNumber,amountData,'b-',label='Original Data')
plt.plot(sampleNumber,fExponential(sampleNumber,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('time(s)')
plt.savefig(fname='query1m.eps',dpi=600,format='eps')
plt.show()

#%%
popt, dump = curve_fit(fExponential,sampleNumber,query1mPerRecord)
print(popt)
plt.plot(sampleNumber,query1mPerRecord,'b-',label='Original Data')
plt.plot(sampleNumber,fExponential(sampleNumber,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('time(s)')
plt.savefig(fname='query1mPerRecord.eps',dpi=600,format='eps')
plt.show()

#%%
popt, dump = curve_fit(fExponential,sampleNumber,queryPerRecord)
print(popt)
plt.plot(sampleNumber,queryPerRecord,'b-',label='Original Data')
plt.plot(sampleNumber,fExponential(sampleNumber,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('average query time(s)')
plt.savefig(fname='queryPerRecord.eps',dpi=600,format='eps')
plt.show()

#%%
popt, dump = curve_fit(fExponential,sampleNumber,query1variantTime)
print(popt)
plt.plot(sampleNumber,query1variantTime,'b-',label='Original Data')
plt.plot(sampleNumber,fExponential(sampleNumber,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('time(s)')
plt.savefig(fname='query1variant.eps',dpi=600,format='eps')
plt.show()

#%%

popt, dump = curve_fit(fNlogN,sampleNumber,query1sampleTime)
print(popt)
plt.plot(sampleNumber,query1sampleTime,'b-',label='Original Data')
plt.plot(sampleNumber,fNlogN(sampleNumber,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('time(s)')
plt.savefig(fname='query1sample.eps',dpi=600,format='eps')
plt.show()

#%%


popt, dump = curve_fit(fExponential,sampleNumber,variantNumber)
print(popt)
plt.plot(sampleNumber,variantNumber,'b-',label='Original Data')
plt.plot(sampleNumber,fExponential(sampleNumber,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('time(s)')
plt.savefig(fname='variantNumber.eps',dpi=600,format='eps')
plt.show()

#%%


#%%
popt, dump = curve_fit(fExponential,sampleNumber,phase1)
print(popt)
plt.plot(sampleNumber,phase1,'b-',label='Original Data')
plt.plot(sampleNumber,fExponential(sampleNumber,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('time(ms)')
plt.savefig(fname='phase1.eps',dpi=600,format='eps')
plt.show()

#%%
popt, dump = curve_fit(fExponential,sampleNumber,phase2)
print(popt)
plt.plot(sampleNumber,phase2,'b-',label='Original Data')
plt.plot(sampleNumber,fExponential(sampleNumber,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('time(ms)')
plt.savefig(fname='phase2.eps',dpi=600,format='eps')
plt.show()

#%%
popt, dump = curve_fit(fExponential,sampleNumber,phase3)
print(popt)
plt.plot(sampleNumber,phase3,'b-',label='Original Data')
plt.plot(sampleNumber,fExponential(sampleNumber,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('time(ms)')
plt.savefig(fname='phase3.eps',dpi=600,format='eps')
plt.show()

#%%
popt, dump = curve_fit(fExponential,amountData,phase3)
print(popt)
plt.plot(amountData,phase3,'b-',label='Original Data')
plt.plot(amountData,fExponential(amountData,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('time(ms)')
plt.savefig(fname='phase3.eps',dpi=600,format='eps')
plt.show()

#%%
popt, dump = curve_fit(fExponential,sampleNumber,phase4)
print(popt)
plt.plot(sampleNumber,phase4,'b-',label='Original Data')
plt.plot(sampleNumber,fExponential(sampleNumber,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('time(ms)')
plt.savefig(fname='phase4.eps',dpi=600,format='eps')
plt.show()

#%%
popt, dump = curve_fit(fExponential,sampleNumber,phase5)
print(popt)
plt.plot(sampleNumber,phase5,'b-',label='Original Data')
plt.plot(sampleNumber,fExponential(sampleNumber,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('time(ms)')
plt.savefig(fname='phase5.eps',dpi=600,format='eps')
plt.show()

#%%


#%%
popt, dump = curve_fit(fExponential,sampleNumber,permutation)
print(popt)
plt.plot(sampleNumber,permutation,'b-',label='Original Data')
plt.plot(sampleNumber,fExponential(sampleNumber,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('time(ms)')
plt.savefig(fname='permutation.eps',dpi=600,format='eps')
plt.show()

#%%
popt, dump = curve_fit(fExponential,sampleNumber,perm1mPerRecord)
print(popt)
plt.plot(sampleNumber,perm1mPerRecord,'b-',label='Original Data')
plt.plot(sampleNumber,fExponential(sampleNumber,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('time(ms)')
plt.savefig(fname='permutation.eps',dpi=600,format='eps')
plt.show()


#%%
popt, dump = curve_fit(fExponential,amountData,permutation)
print(popt)
plt.plot(amountData,permutation,'b-',label='Original Data')
plt.plot(amountData,fExponential(amountData,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('amount of data')
plt.ylabel('time(ms)')
plt.savefig(fname='permutation.eps',dpi=600,format='eps')
plt.show()

#%%
popt, dump = curve_fit(fExponential,sampleNumber,real)
print(popt)
plt.plot(sampleNumber,real,'b-',label='Original Data')
plt.plot(sampleNumber,fExponential(sampleNumber,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('time(ms)')
plt.savefig(fname='real.eps',dpi=600,format='eps')
plt.show()

#%%
popt, dump = curve_fit(fExponential,sampleNumber, compresSize)
print(popt)
plt.plot(sampleNumber,compresSize,'b-',label='Original Data')
plt.plot(sampleNumber,fExponential(sampleNumber,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('compressed size')
plt.savefig(fname='compresSize.eps',dpi=600,format='eps')
plt.show()

#%%
popt, dump = curve_fit(fExponential,sampleNumber,iteration)
print(popt)
plt.plot(sampleNumber,iteration,'b-',label='Original Data')
plt.plot(sampleNumber,fExponential(sampleNumber,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('times')
plt.savefig(fname='iteration.eps',dpi=600,format='eps')
plt.show()

#%%
popt, dump = curve_fit(fExponential,sampleNumber,permPerIter)
print(popt)
plt.plot(sampleNumber,permPerIter,'b-',label='Original Data')
plt.plot(sampleNumber,fExponential(sampleNumber,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('time(ms)')
plt.savefig(fname='permPerIter.eps',dpi=600,format='eps')
plt.show()

#%%
popt, dump = curve_fit(fExponential,sampleNumber[0:bitVectorNumber.shape[0]],bitVectorNumber)
print(popt)
plt.plot(sampleNumber[0:bitVectorNumber.shape[0]],bitVectorNumber,'b-',label='Original Data')
plt.plot(sampleNumber[0:bitVectorNumber.shape[0]],fExponential(sampleNumber[0:bitVectorNumber.shape[0]],*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('number of bit vectors')
plt.savefig(fname='bitVectorNumber.eps',dpi=600,format='eps')
plt.show()

#%%
variantNumber = variantNumber[0:bitVectorNumber.shape[0]]
popt, dump = curve_fit(fExponential,variantNumber,bitVectorNumber)
print(popt)
plt.plot(variantNumber,bitVectorNumber,'b-',label='Original Data')
plt.plot(variantNumber,fExponential(variantNumber,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of records')
plt.ylabel('number of bit vectors')
plt.savefig(fname='dump.eps',dpi=600,format='eps')
plt.show()

#%%
popt, dump = curve_fit(fExponential,sampleNumber,rest1mPerRecord)
print(popt)
plt.plot(sampleNumber,rest1mPerRecord,'b-',label='Original Data')
plt.plot(sampleNumber,fExponential(sampleNumber,*popt),'r-',label='Fitted Curve')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('number of bit vectors')
plt.savefig(fname='bitVectorNumber.eps',dpi=600,format='eps')
plt.show()

#%%
plt.plot(amountData,perm1m,'b-',label='Permutation time of original GTC')
plt.plot(amountData,perm1mNew,'r-',label='Permutation time of changed GTC')
plt.legend()
plt.xlabel('amount of data')
plt.ylabel('time(s)')
plt.savefig(fname='permImprovment.eps',dpi=600,format='eps')
plt.show()

#%%
plt.plot(sampleNumber,real1m,'b-',label='Real time of original GTC')
plt.plot(sampleNumber,real1mNew,'r-',label='Real time of changed GTC')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('time(s)')
plt.savefig(fname='realImprovment.eps',dpi=600,format='eps')
plt.show()

#%%
plt.plot(sampleNumber,real1mPerRecord,'b-',label='original GTC')
plt.plot(sampleNumber,real1mNewPerRecord,'r-',label='changed GTC')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('average query time for one record(s)')
plt.savefig(fname='realImprovment.eps',dpi=600,format='eps')
plt.show()

#%%
plt.plot(sampleNumber,perm1mPerRecord,'b-',label='original GTC')
plt.plot(sampleNumber,perm1mNewPerRecord,'r-',label='changed GTC')
plt.legend()
plt.xlabel('number of samples')
plt.ylabel('average query time for one record(s)')
plt.savefig(fname='realImprovment.eps',dpi=600,format='eps')
plt.show()

#%%
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
percentage = perm1mPerRecord/real1mPerRecord
plt.plot(sampleNumber,percentage,'b-')
plt.xlabel('number of samples')
plt.ylabel('average query time for one record(s)')
plt.savefig(fname='percentage.eps',dpi=600,format='eps')
plt.show()