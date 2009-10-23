from defs_moduleByTableLookUp import *
from JamesDefs import groupByField
from doug_hypergeometric import hyperGeoPvalue
from time import time
from sets import Set




#========================= User Defined Variables =========================
inputMotifComboList = '/Users/biggus/Documents/MBGB/Rotations/James/Data/motifCombos/tissue_2Combos.0_1culled.nr.txt'
originalMotifs      = '/Users/biggus/Documents/MBGB/Rotations/James/Data/DegenMotifs/aedesAnopheles7mer.TSS.nr.motif'
matrixFile          = '/Users/biggus/Documents/MBGB/Rotations/James/Data/MotifMaps/motif_2Combos.txt'

outFile = ''
#==========================================================================

inputMotifComboList = map(lambda line: line.strip().split('\t'), open(inputMotifComboList, 'rU').readlines())

originalMotifs      = map(lambda line: line.strip(), open(originalMotifs, 'rU').readlines())

matrixFile          = readInMatrixFromFile(matrixFile)

setOfLargerCombos = []



for combo in inputMotifComboList:
    
    for motif in originalMotifs:
        
        if motif not in combo:
            
            setOfLargerCombos.append(Set(combo+[motif]))
            
            
            
#  convert to set and back to ensure no repetition

lenSet = len(setOfLargerCombos)
setOfLargerCombos = Set(setOfLargerCombos)
lenSet = len(setOfLargerCombos)

x=0