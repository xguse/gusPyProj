from sets import Set
from time import time

#========================= User Defined Variables =========================
pValsFile      = '/Users/biggus/Documents/MBGB/Rotations/James/Data/2KB/2kb_ClusterHyperGeo/comboHyperGeo/12_time-course*3h_UP.BH1.pvals.txt'

#===============================================

pValsFile = map(lambda line: line.strip(), open(pValsFile, 'rU').readlines())

nrListOfMotifs = []

for each in pValsFile:
    motifs = each.split('\t')[0].split('_')
    nrListOfMotifs+= motifs
    
nrListOfMotifs = list(Set(nrListOfMotifs))
nrListOfMotifs.sort()

for i in nrListOfMotifs:
    print i