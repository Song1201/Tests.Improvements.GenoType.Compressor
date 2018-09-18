#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 15:19:44 2018

@author: songchen
"""

import vcf
import BitVector

#%%
def bitVectorsToVcf(bitVectors,POSs,refs,alts,template,outputFile):
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
            record.samples[sampleNumber].gt_nums = ''
            
    