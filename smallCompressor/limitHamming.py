#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 25 22:47:20 2018
Limited Hamming
@author: songchen
"""

#%%
import vcf
import BitVector as bv
import numpy as np
import vcfpy as vp
import compressLib as cl
import time
import sys

startPos = 16800001
endPos = 17000000

bitVectors,POSs,refs,alts,sampleNames = cl.readVcfgzToBitVectors(filename
                                        ='gene1000chr22.vcf.gz',startPos=startPos,
                                        endPos=endPos)

colBitVecs = cl.columnToRow(bitVectors)
#%%
# Inputs bitVectors
bitVectors = colBitVecs

# Record Arrays
allUnfitPos = np.array([])
allFitPos = np.array([])
allMoves = np.array([])
allTotalHams = np.array([])


startTime = time.time()

# Setup basic variables
amountBitVecs = len(bitVectors)
length = bitVectors[0].length()
limit = int(0.2*amountBitVecs) # Move at most 0.2 of all vectors
# The initial searching depth. After searching depth*limit without finding a good switch, the searching is reset and the threshold will be lower.
depth = 0.2
finalDepth = 1
depthGrowth = 0.1# How much the depth of searching grows, once the searching is reset
# An imperical setting, threshold for the profit rate, profit/cost
iniThresh = int(length/600)*50 
# How much the threshold decreases, once the searching is reset
threshDecay = int(depthGrowth/2*iniThresh) 
IDs = list(range(amountBitVecs)) # Initialize IDs
amountLevel = int((finalDepth-depth)/depthGrowth+1) # Amount of sleep levels
print(iniThresh, limit)

# Prepare for searching
moves = 0
# current threshold is the threshold hold used for this sleep level, which means it won't be changed if the fail counter is not reset. 
currThresh = iniThresh
activePos = {}
sleepPos = [{} for i in range(amountLevel)]
sleepLevel = 0
# Setup the lowest threshold, so that once the profit rate is lower than this value, put it to the last sleep level.
lowestThresh = iniThresh - threshDecay*(amountLevel-1)
# Use a numpy array to store all the hamming distances between 2 bit vectors
allHamsNext = np.array([bitVectors[ID].hamming_distance(bitVectors[ID+1]) for ID in
                  IDs[0:-1]],dtype=np.int16)

totalHam = np.sum(allHamsNext)
# At the beginning all IDs are movable, when a bit vector is marked as unmovalbe, set the corresponding value to -1, so that it won't be again selected as the unfit vector.


# when the amount of moved bit vectors are still under the limit
while moves < limit:
    allHamsBoth = np.append(allHamsNext,0)
    allHamsBoth[1:] += allHamsNext
    # Rank all the positions according to their allHamsBoth, in an descending order. When search for a good switch, search from the left, assuming the unfitting one usually has a large Hamming distance from both sides.
    POSs = np.flip(allHamsBoth.argsort(),-1)
    # threshold will change during the switch searching for the same position, updating to the best profit rate, in order to find the best switch and be set back to the currThresh once a searching for the best switch for a position is finished.
    threshold = currThresh
    # FailCounter counter for how many times a switch can't be found, also use to mark where the search should begin
    failCounter = 0
    for unfitPos in POSs:
        # Add active postions
        if unfitPos not in activePos:
            activePos[unfitPos] = {a for a in POSs[failCounter+1:]}  
        # Define the searching area, don't have to be copy, use copy to avoid potential problem
        searchArea = activePos[unfitPos].copy()
        found = False # A flag for whether find a good switch or not

        updateHams = []
        updatePos = [] # Use list because set has an order chaos
        for pos in searchArea:
            # Assign the left and right position
            if pos < unfitPos:
                leftPos,rightPos = pos,unfitPos
            else:
                leftPos,rightPos = unfitPos,pos 
            
            
            # The if else branches below are used to calculate the whose Hamming distance to the next position need to be updated and update to what. The profit will be calculated at last.
            # If they are besides each other
            if rightPos - leftPos == 1:
                tempUpdatePos = {i for i in [leftPos-1,rightPos] if -1<i<amountBitVecs-1}
                if leftPos == 0: 
                    tempUpdateHams = [(bitVectors[IDs[0]].hamming_distance(
                            bitVectors[IDs[2]]))]                
                    #profit = allHamsNext[1] - tempUpdateHams[0]
                elif rightPos == amountBitVecs-1:
                    tempUpdateHams=[bitVectors[IDs[-3]].hamming_distance(bitVectors[IDs[-1]])]
                    #profit = allHamsNext[-2] - tempHams[0]
                else:
                    tempUpdateHams = [bitVectors[IDs[leftPos-1]].hamming_distance(
                            bitVectors[IDs[rightPos]])] # left and right switch
                    tempUpdateHams.append(bitVectors[IDs[leftPos]].hamming_distance(
                            bitVectors[IDs[rightPos+1]]))
                    #profit = allHamsNext[leftPos-1] + allHamsNext[rightPos] - sum(tempHams)

            # If they are separated by 1 position
            elif rightPos - leftPos == 2:
                tempUpdatePos = {i for i in list(range(leftPos-1,rightPos+1)) if                          
                                 -1<i<amountBitVecs-1}
                if leftPos == 0:
                    tempUpdateHams = [allHamsNext[1]]
                    tempUpdateHams.append(allHamsNext[0])
                    tempUpdateHams.append(bitVectors[IDs[0]].hamming_distance(
                            bitVectors[IDs[3]]))
                    #profit = allHamsNext[2] - tempHams[-1]
                elif rightPos == amountBitVecs-1:
                    tempUpdateHams = [bitVectors[IDs[-4]].hamming_distance(
                            bitVectors[IDs[-1]])]
                    tempUpdateHams.append(allHamsNext[-1])
                    tempUpdateHams.append(allHamsNext[-2])
                    #profit = allHamsNext[-3] - tempHams[0]
                else:
                    tempUpdateHams = [bitVectors[IDs[leftPos-1]].hamming_distance(
                            bitVectors[IDs[rightPos]])]
                    tempUpdateHams.append(allHamsNext[leftPos+1])
                    tempUpdateHams.append(allHamsNext[leftPos])
                    tempUpdateHams.append(bitVectors[IDs[leftPos]].hamming_distance(
                            bitVectors[IDs[rightPos+1]]))
                    #profit = allHamsNext[leftPos-1] + allHamsNext[rightPos] - tempHams[0] - tempHams[-1]
            # If they are separated by more than 1 position
            else: 
                tempUpdatePos = {i for i in [leftPos-1,leftPos,rightPos-1,rightPos] if 
                                 -1<i<amountBitVecs-1}
                if leftPos == 0:
                    tempUpdateHams = [bitVectors[IDs[rightPos]].hamming_distance(
                            bitVectors[IDs[1]])]
                else:
                    tempUpdateHams = [bitVectors[IDs[leftPos-1]].hamming_distance(
                            bitVectors[IDs[rightPos]])]
                    tempUpdateHams.append(bitVectors[IDs[rightPos]].hamming_distance(
                            bitVectors[IDs[leftPos+1]]))
                if rightPos == amountBitVecs-1:
                    tempUpdateHams.append(bitVectors[IDs[-2]].hamming_distance(
                            bitVectors[IDs[leftPos]]))
                else:
                    tempUpdateHams.append(bitVectors[IDs[rightPos-1]].hamming_distance(
                            bitVectors[IDs[leftPos]]))
                    tempUpdateHams.append(bitVectors[IDs[leftPos]].hamming_distance(
                            bitVectors[IDs[rightPos+1]]))
                    
            profit = np.sum(allHamsNext[list(tempUpdatePos)] - tempUpdateHams)
            # Deside switch or not according to profit and cost    
            if profit < 0:
                activePos[unfitPos].remove(pos)
                continue
            cost = int(IDs[unfitPos]==unfitPos) + int(IDs[pos]==pos) - int(IDs[
                    unfitPos]==pos) - int(IDs[pos]==unfitPos)   
            if cost < 0:
                # Select this switch as the best switch and do it anyway
                updateHams = tempUpdateHams
                updatePos = list(tempUpdatePos)
                fitPos = pos
                costSwitch = cost
                found = True
                #whySwitch = 'Minus Cost'
                break
            elif profit==0:
                activePos[unfitPos].remove(pos)
                continue
            elif cost==0: # profit > 0, no cost, do the switch
                updateHams = tempUpdateHams
                updatePos = list(tempUpdatePos)
                fitPos = pos
                costSwitch = cost
                found = True
                #whySwitch = 'No cost'
                break    
            else:
                rate = profit/cost
            if rate > threshold: 
                updateHams = tempUpdateHams
                updatePos = list(tempUpdatePos)
                fitPos = pos
                threshold = rate
                costSwitch = cost
                found = True
                continue
            # The switch is even worse than the current threshold, put the switch to sleep, the sleep level is based on its profit/cost
            elif lowestThresh < rate <= currThresh: 
                level = int((iniThresh-rate)/threshDecay)
                activePos[unfitPos].remove(pos)
                if unfitPos in sleepPos[level]:
                    sleepPos[level][unfitPos].add(pos)
                else:
                    sleepPos[level][unfitPos] = {pos}
                continue
            elif rate <= lowestThresh: # Rate is too low, put to the deepest sleep area
                activePos[unfitPos].remove(pos)
                if unfitPos in sleepPos[-1]:
                    sleepPos[-1][unfitPos].add(pos)
                else:
                    sleepPos[-1][unfitPos] = {pos}
                continue
            else: # The switch is not the best, but good enough to be active
                continue
                    
        
        if found:
            #Switch                
            allUnfitPos = np.append(allUnfitPos,unfitPos)
            allFitPos = np.append(allFitPos,fitPos)
            IDs[fitPos],IDs[unfitPos] = IDs[unfitPos],IDs[fitPos]
            moves += costSwitch
            allMoves = np.append(allMoves,moves)
            updatePos.sort() # Because when convert from set, the order maybe weird
            allHamsNext[updatePos] = updateHams
#            for i in range(amountBitVecs-1):
#                if allHamsNext[i] != bitVectors[IDs[i]].hamming_distance(bitVectors[IDs[i+1]]):
#                    print('Hamming Wrong!',i,unfitPos,fitPos)
#                    sys.exit()

            
            
            # the positions that can be reactivated
            refreshArea = {i for i in list(range(unfitPos-1,unfitPos+2))+list(
                    range(fitPos-1,fitPos+2)) if -1 < i < amountBitVecs}
            if abs(unfitPos - fitPos) == 2:
                refreshArea.remove((unfitPos+fitPos)/2)
            
            # Update active positions
            for key in activePos:
                activePos[key].update(refreshArea)
                if key in refreshArea:
                    activePos[key].remove(key)

            newTotalHam = np.sum(allHamsNext)
            allTotalHams = np.append(allTotalHams,newTotalHam)
            deltaHam = totalHam - newTotalHam
            totalHam = newTotalHam
            
            print('Total Hamming Distance:',totalHam,' Moves:',moves,' profit:',deltaHam,' FailCounter:',failCounter,' Current Thresh:',currThresh)
            break


        else:
            failCounter += 1
            if failCounter % (0.05*amountBitVecs) == 0:
                print('Still Running, failCounter:',failCounter)
            if failCounter > depth*limit:
                failCounter = 0
                currThresh -= threshDecay
                depth += depthGrowth
                if depth >= 1: #  Searching is too difficult, just finish ASAP
                    currThresh = 0
                    continue
                # Wake up the sleeping positions in this level
                for key in sleepPos[sleepLevel]:
                    activePos[key].update(sleepPos[sleepLevel][key])
                sleepPos[sleepLevel].clear()
                sleepLevel += 1
                print('Bar too high, reset threshold. Sleep Level:',sleepLevel)
                break
                
endTime = time.time()
print('Runtime:',endTime-startTime)
            
        




allUnfitPos.dump('unfitPos'+str(startPos)[:4]+'_'+str(endPos)[0:4]+'_50.dat')
allFitPos.dump('fitPos'+str(startPos)[:4]+'_'+str(endPos)[0:4]+'_50.dat')
allMoves.dump('moves'+str(startPos)[:4]+'_'+str(endPos)[0:4]+'_50.dat')
allTotalHams.dump('totalHams'+str(startPos)[:4]+'_'+str(endPos)[0:4]+'_50.dat')
