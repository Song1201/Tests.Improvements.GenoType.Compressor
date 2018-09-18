#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 18:44:00 2018

@author: songchen
"""

import vcf
import numpy as np
import BitVector as bv
import vcfpy as vp
import compressLib as cl
import importlib 
importlib.reload(cl)
import matplotlib.pyplot as plt
import os 

#%% Get the first record in PyVCF in simulated population 1
readerPyVcf = vcf.Reader(filename='simulatedPopulation1.vcf.gz')
recordPyVcf = next(readerPyVcf)
print(recordPyVcf)
print(recordPyVcf.get_hets())

#%% Get the first record in VcfPy in simulated population 1
readerVcfPy = vp.Reader.from_path('simulatedPopulation1.vcf.gz')
recordVcfPy = next(readerVcfPy)
a = recordVcfPy.calls[278].gt_alleles
print(a)

#%% Check the template
reader = vcf.Reader(filename='template.vcf')
for record in reader:
    print(record)
    print(record.get_hets())

#%% Get the first severl records in simulatedPopulation and another file
several = 5    
readerOrigin = vcf.Reader(filename='simulatedPopulation1_2.vcf.gz')
readerChange = vcf.Reader(filename='changed.vcf')
for i in range(several):
    recordOrigin = next(readerOrigin)
    print(recordOrigin)  
    print(recordOrigin.get_hets())
    recordChange = next(readerChange)
    print(recordChange)
    print(recordChange.get_hets())
    

#%%
reader = vp.Reader.from_path('template.vcf')
writer = vp.Writer.from_path('tryOutput.vcf',reader.header)
record = next(reader)
sub2 = vp.Substitution(type_='SNV',value='Z')
record.ALT.append(sub2)
#call1 = record.calls[278]
record.call_for_sample = call1
#record.calls[278].gt_alleles.clear()
#record.calls[278].gt_alleles.append(2)
#record.calls[278].gt_alleles.append(0)
print(record.calls[278].gt_alleles)
writer.write_record(record)
writer.close()

#%% Check the result in tryOutput.vcf
readerTryOutput = vcf.Reader(filename='tryOutput.vcf')
recordTryOutput = next(readerTryOutput)
print(recordTryOutput)
print(recordTryOutput.get_hets())

#%% Test for bitVectorsToVcf
# inputs
template = 'template.vcf'

reader = vcf.Reader(filename=template)
record = next(reader)
writer = vcf.Writer(open(outputFile,'w'),reader)
for recordNumber in range(len(POSs)):
    record.POS = POSs[recordNumber]
    record.REF = refs[recordNumber]
    record.ALT = alts[recordNumber]
    for bitVectorNumber in range(len(record.ALT)-1,-1,-1):
        bitVectors[bitVectorNumber]
    genotypes = bitVectors
    for sampleNumber in range(len(record.samples)):
        record.samples[sampleNumber].gt_nums = 'dd'
        
#%% Test for saveTemplate
# Inputs
filename = 'simulatedPopulation1.vcf.gz'
template = 'template.vcf'        
        
# Use pyVCF to get the samples whose genotypes are not 0/0
reader = vcf.Reader(filename=filename)
record = next(reader)
hets = record.get_hets()
hetSamples = []
for call in hets:
    hetSamples.append(call.sample)
# Use VCFpy to write the record to all 0/0,clear the ALT, set the REF to empty string,CHROM to  0 and save the template    
reader = vp.Reader.from_path(path=filename)
record = next(reader)
record.CHROM = 0
record.POS=0
record.REF=''
record.ALT.clear()
for sample in hetSamples:
    record.call_for_sample[sample] = vp.Call(sample=sample,data={'GT': '0/0'})
writer = vp.Writer.from_path(path=template,header=reader.header)
writer.write_record(record)
writer.close()

#%% Test for readVcfgzToBitVectors
# Input: filename, startPos, endPos,chrom
filename = 'gene1000chr22.vcf.gz'
startPos = 16050075
endPos = 16500000
chrom = '22'

reader = vcf.Reader(filename=filename)
pop = len(reader.samples) # population size, the number of samples
bitVectors = [] # a list for restoring the bit vectors, don't use np.array because restore bit vectors into numpy array will convert the bit vectors to floats, which may lead to huge memory consume.
sampleNames = reader.samples
# using lists because list use smaller memory, though sacrifised the speed.
POSs = []
refs = []
alts = [] 
for record in reader.fetch(chrom,startPos,endPos): 
    #ALT read from pyVcf is a list of object, so need to do some process
    if len(record.ALT) == 1:
        alts.append(record.ALT[0].sequence)
        POSs.append(record.POS)
        refs.append(record.REF)
    else:
        for alt in record.ALT:
            if alt.type == 'SNV':
                POSs.append(record.POS)
                refs.append(record.REF)
                alts.append(alt.sequence)
            else:
                POSs.append(record.POS)
                refs.append(record.REF)
                alts.append(alt.type)
    # refer the GTC paper to see why the length of bit vectors equal to the twice of population size
    recordBitVectors = [] # bit vectors for this record
    for i in range(len(record.ALT)):
        recordBitVectors.append(bv.BitVector(size=2*pop))
    # can't handle the situation for unknown variants.
    for sample in record.get_hets():
        sampleNumber = sampleNames.index(sample.sample)
        if sample.gt_nums[0] != '0':
            recordBitVectors[int(sample.gt_nums[0])-1][sampleNumber*2] = 1
        if sample.gt_nums[2] != '0':
            recordBitVectors[int(sample.gt_nums[2])-1][sampleNumber*2+1] = 1
    for recordBitVector in recordBitVectors:
        bitVectors.append(recordBitVector)

#%% Test for columnToRow
# Inputs: bitVectors


colBitVecs = []
height = len(bitVectors) # the height of the input bit array
width = bitVectors[0].length()
# Build a list of bit vectors, each bit vector is a column in the iput bit array
for colNumber in range(width):
    colBitVecs.append(bv.BitVector(size=height))
# Search for 1s in the rows of input bit array, and set the corresponding bits in the new bit vectors    
for rowNumber in range(height):
    setPos = 0 # Initialize the position of the set bit
    # When it comes to no more set bit in the bit vector, the method .next_set_bit() returns -1
    while True : 
        setPos = bitVectors[rowNumber].next_set_bit(setPos)
        if setPos == -1:
            break
        colBitVecs[setPos][rowNumber] = 1
        setPos += 1
    
#%% To see if there are unknown genotypes in a file.
reader = vcf.Reader(filename='gene1000chr22.vcf.gz')
i = 0
for record in reader:
    i += 1
    if i % 10000 == 0:
        print(i)
    if record.get_unknowns() == []:
        ;
    else:
        print(record)
        break
    
#%% Test for heurisHamming
# Inputs: bitVectors
bitVectors = colBitVecs

positions = range(len(bitVectors)) # Original position
# The final ordered bit vectors, use the first bit vector as the default start position
orderedBitVecs = [bitVectors[0].deep_copy()] # 
allHamming = 0 # The sum of all hamming distance
# Initialize the rest positions
restPos = list(positions)
restPos.remove(0)
# Initialize the nearest position to the start position
nearestPos = 0 
# loop when there are still unsearched positions
while len(restPos)>0:          
    # Use the longest possible Hamming distance to initialize
    hamming = bitVectors[0].length() 
    for pos in restPos:
        newHamming = orderedBitVecs[-1].hamming_distance(bitVectors[pos])
        if newHamming < hamming:
            hamming = newHamming
            nearestPos = pos
    allHamming += hamming 
    orderedBitVecs.append(bitVectors[nearestPos].deep_copy())
    restPos.remove(nearestPos)

    print(allHamming,'',len(orderedBitVecs))
        
#%% Test of limitHamming
# Inputs bitVectors
bitVectors = colBitVecs
limit = 300
iniThresh = 100

# Outputs

threshold = iniThresh
IDs = list(range(len(bitVectors))) # Initialize IDs
positions = list(range(len(bitVectors))) # Just a list from 0 to the last
#IDs = np.arange(len(bitVectors),dtype=np.int16) # Initialize IDs
amountBitVecs = len(bitVectors)
length = bitVectors[0].length()
# Use a numpy array to store all the hamming distances from both side (or one side for the first and the last bit vector) for all positions
allHamsBoth = np.zeros(amountBitVecs,dtype=np.int16)
allHamsNext = np.zeros(amountBitVecs-1,dtype=np.int16) 
# In this loop, IDs are in the original order, which means IDs[i]==i
for ID in IDs:
    if ID == 0:
        allHamsNext[ID]=allHamsBoth[ID] = bitVectors[ID].hamming_distance(bitVectors[ID+1])        
    elif ID == amountBitVecs-1:
        allHamsBoth[ID] = allHamsNext[-1]
    else:
        allHamsNext[ID] = bitVectors[ID].hamming_distance(bitVectors[ID+1])
        allHamsBoth[ID] = allHamsNext[ID-1] + allHamsNext[ID]

movableHams = allHamsBoth # At the beginning all IDs are movable, when a bit vector is marked as unmovalbe, set the corresponding value to -1, so that it won't be again selected as the unfit vector.
failCounter = 0
# when the amount of moved bit vectors are still under the limit
while limit > 0:
    # unfitPos will never be the last position
    unfitPos = movableHams.argmax()
    fitPos = unfitPos # Initialize the fitPos to the impossible one, for later verify
    # Store the value and the index which need to be updated in the allHamsNext,allHamsBoth
    betterHams = []
    betterPos = []
    effectBoth = [] # The value will be calculated at last
    # The positions in which we will search for a good switch
    positions.remove(unfitPos)
    for pos in positions:
        # Temporary storage, if it is good enough, update them to value and index
        tempHams = []
        tempPos = []
        tempBoth = []
        # Assign the left and right position
        if pos < unfitPos:
            leftPos = pos
            rightPos = unfitPos
        else:
            leftPos = unfitPos
            rightPos = pos
        # If they are besides each other
        if rightPos - leftPos == 1:
            if leftPos == 0: 
                tempPos.append(1)
                tempHams.append(bitVectors[IDs[0]].hamming_distance(bitVectors[IDs[2]]))                
                profit = allHamsNext[1] - tempHams[0]
                tempBoth = [1,2]
            elif rightPos == amountBitVecs-1:
                tempPos.append(-2)
                tempHams.append(bitVectors[IDs[-3]].hamming_distance(bitVectors[IDs[-1]]))                
                profit = allHamsNext[-3] - tempHams[0]
                tempBoth = [-3,-2]
            else:
                tempHams.append(bitVectors[IDs[leftPos-1]].hamming_distance(
                        bitVectors[IDs[rightPos]])) # left and right switch
                tempHams.append(bitVectors[IDs[leftPos]].hamming_distance(
                        bitVectors[IDs[rightPos+1]]))
                tempPos.append(leftPos-1)
                tempPos.append(rightPos)
                profit = allHamsNext[leftPos-1] + allHamsNext[rightPos] - sum(tempHams)
                tempBoth = [leftPos-1,leftPos,rightPos,rightPos+1]
        # If they are separated by 1 position
        elif rightPos - leftPos == 2:
            if leftPos == 0:
                tempPos = [0,1,2]
                tempHams.append(allHamsNext[1])
                tempHams.append(allHamsNext[0])
                tempHams.append(bitVectors[IDs[0]].hamming_distance(bitVectors[IDs[3]]))
                profit = allHamsNext[2] - tempHams[-1]
                tempBoth = [0,2,3]
            elif rightPos == amountBitVecs-1:
                tempPos = [-3,-2,-1]
                tempHams.append(bitVectors[IDs[-4]].hamming_distance(bitVectors[IDs[-1]]))
                tempHams.append(allHamsNext[-1])
                tempHams.append(allHamsNext[-2])
                profit = allHamsNext[-3] - tempHams[0]
                tempBoth = [-4,-3,-1]
            else:
                tempPos.append(leftPos-1)
                tempHams.append(bitVectors[IDs[leftPos-1]].hamming_distance(
                        bitVectors[IDs[rightPos]]))
                tempPos.append(leftPos)
                tempHams.append(allHamsNext[leftPos+1])
                tempPos.append(leftPos+1)
                tempHams.append(allHamsNext[leftPos])
                tempPos.append(rightPos)
                tempHams.append(bitVectors[IDs[leftPos]].hamming_distance(
                        bitVectors[IDs[rightPos+1]]))
                profit = allHamsNext[leftPos-1] + allHamsNext[rightPos] - tempHams[0] - tempHams[-1]
                tempBoth = [leftPos-1,leftPos,rightPos,rightPos+1]
        else: 
            if leftPos == 0:
                tempPos.append(0)
                tempHams.append(bitVectors[IDs[rightPos]].hamming_distance(
                        bitVectors[IDs[1]]))
                tempBoth = [0,1]
            else:
                tempPos.append(leftPos-1)
                tempHams.append(bitVectors[IDs[leftPos-1]].hamming_distance(
                        bitVectors[IDs[rightPos]]))
                tempPos.append(leftPos)
                tempHams.append(bitVectors[IDs[rightPos]].hamming_distance(
                        bitVectors[IDs[leftPos+1]]))
                tempBoth = [leftPos-1,leftPos,leftPos+1]
            if rightPos == amountBitVecs-1:
                tempPos.append(-1)
                tempHams.append(bitVectors[IDs[-2]].hamming_distance(
                        bitVectors[IDs[leftPos]]))
                tempBoth.append(-2)
                tempBoth.append(-1)
            else:
                tempPos.append(rightPos-1)
                tempHams.append(bitVectors[IDs[rightPos-1]].hamming_distance(
                        bitVectors[IDs[leftPos]]))
                tempPos.append(rightPos)
                tempHams.append(bitVectors[IDs[leftPos]].hamming_distance(
                        bitVectors[IDs[rightPos+1]]))
                tempBoth.append(rightPos-1)
                tempBoth.append(rightPos)
                tempBoth.append(rightPos+1)
            profit = np.sum(allHamsNext[tempPos]) - sum(tempHams)
            
        if profit < 0:
            continue
        cost = int(IDs[leftPos]==leftPos) + int(IDs[rightPos]==rightPos) - int(IDs[
                leftPos]==rightPos) - int(IDs[rightPos]==leftPos)
        if cost < 0:
            # Select this switch as the better switch and do it anyway
            betterHams = tempHams
            betterPos = tempPos
            effectBoth = tempBoth
            fitPos = pos
            break
        elif profit==0:
            continue
        elif cost==0: # profit > 0, no cost, do the switch
            betterHams = tempHams
            betterPos = tempPos
            effectBoth = tempBoth
            fitPos = pos
            break
        elif profit/cost > threshold: 
            betterHams = tempHams
            betterPos = tempPos
            fitPos = pos
            threshold = profit/cost
            continue
        else:
            continue
    if fitPos == unfitPos: # No fit position has been found, look for the next unfit position
        movableHams[unfitPos] = -1 # So that it won't be located again
        failCounter += 1
#        if failCounter % 20 == 0:
#            threshold -= 2
        print('No fitting position found.')
    else: # A fitting position has been found
        # Do the switch
        dump = IDs[fitPos]
        IDs[fitPos] = IDs[unfitPos]
        IDs[unfitPos] = dump
        limit -= cost
        allHamsNext[betterPos] = betterHams
        failCounter = 0
        threshold = iniThresh
        for pos in effectBoth:
            if pos==0:
                allHamsBoth[pos] = allHamsNext[0]
            elif pos == amountBitVecs-1 or pos == -1:
                allHamsBoth[pos] = allHamsNext[-1]
            elif pos < 0:
                allHamsBoth[pos] = allHamsNext[pos] + allHamsNext[pos+1]
            else:
                allHamsBoth[pos] = allHamsNext[pos-1] + allHamsNext[pos]
        movableHams = allHamsBoth # Again, all positions can be moved
        positions = list(range(amountBitVecs))
        
        outputVecs = []
        for ID in IDs:
            outputVecs.append(bitVectors[ID].deep_copy())
        newTotalHam = 0
        for i in range(len(outputVecs)-1):
            newTotalHam += outputVecs[i].hamming_distance(outputVecs[i+1])
        totalHam = np.sum(allHamsNext)
        print('Total Hamming Distance:',totalHam,' New Hamming:',newTotalHam,
              '  Rest Move:',limit)
    
                    
#%% test for writeBitVecs
# Input
#bitVectors
i = 0
for vec in bitVectors:
    i += 1
    file = open((os.getcwd() + '/bitVectors/bitVector{0:.0f}.bits'.format(i)),'wb')
    vec.write_to_file(file)
    file.close()
    
#%% Read bit 
rBitVectors = []
for i in range(1,7222):
    a = bv.BitVector(filename=(os.getcwd() +'/bitVectors/bitVector{0:.0f}.bits').format(i))
    rBitVectors.append(a.read_bits_from_file(20000))
    if i % 1000 == 0:
        print(i)

#%%
for i in range(1,7222):
    if i % 1000 == 0:
        print('Checked:',i)
    if bitVectors[i] == rBitVectors[i]:
        pass
    else: 
        print('Wrong! ',i)
        sys.exit()
    
    

#%% Check the statistic characteristic of the Hamming distances
bitVectors = heurisBitVecs
hammings = np.zeros(len(bitVectors)-1,dtype=np.uint16)
for bitVecNumber in range(len(bitVectors)-1):
    hammings[bitVecNumber] = bitVectors[bitVecNumber].hamming_distance(
            bitVectors[bitVecNumber+1])
    
#%% Test the speed of list and dictionary
import time
dictionary = {}
size = 100000
for i in range(size):
    dictionary['S_{0:.0f}'.format(i)] = i

name = 'S_20000'
startTime = time.time()
i = dictionary[name]
endTime = time.time()
print(endTime - startTime)

aList = []
for i in range(size):
    aList.append('S_{0:.0f}'.format(i))

name = 'S_20000'
startTime = time.time()
i = aList.index(name)
endTime = time.time()
print(endTime - startTime)

#%% Test for countNonRef
import BitVector as bv
import vcf
# Input:
filename = 'gene1000chr22.vcf.gz'
startPos = 16050075
endPos = 16500000
chrom = '22'

# Function
reader = vcf.Reader(filename=filename)
pop = len(reader.samples) # population size, the number of samples
bitVectors = [] # a list for restoring the bit vectors, don't use np.array because restore bit vectors into numpy array will convert the bit vectors to floats, which may lead to huge memory consume.
sampleNumber = 0
sampleDict = {}
for sampleName in reader.samples:
    sampleDict[sampleName] = sampleNumber
    sampleNumber += 1

startPos -= 1 # the reader.fetch method doesn't include startPos itself, which is weird
recordNo = 0
for record in reader.fetch(chrom,startPos,endPos):     
    # can't handle the situation for unknown variants.
    bitVec = bv.BitVector(size=pop)
    for call in record.get_hets():
        bitVec[sampleDict[call.sample]] = 1 
    bitVectors.append(bitVec)
    recordNo += 1
    if recordNo % 1000 == 0:
        print('Number of read records:',recordNo)
        
#%% Test for 