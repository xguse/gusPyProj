import string
from gusPyCode.defs import JamesDefs
import sys



#========================= User Defined Variables =========================

# File paths:

# IN:
codingBoundsList       = open('/Users/biggus/Documents/James/Data/Tu_miRNA/Cq_500afterCoding.newCoords.txt','r').readlines()
resolvedConflictsList  = open('/Users/biggus/Documents/James/Data/Tu_miRNA/Cq_500afterCoding.newCoords.fjoinResolved.txt','r').readlines()

# OUT:
reCombinedCoordsFILE   = open('/Users/biggus/Documents/James/Data/Tu_miRNA/Cq_500afterCoding.newCoords.usuable.txt','w')

boundaryRegion = 'down' # 'up' or 'down'   <-<-<-<-<-<-!!!!!!!!!__DO NOT FORGET TO CHANGE THIS__!!!!!!!!!!!


shortestUsableBdryReg  = 15


#==========================================================================

#--------- Script Specific Function Definitions ---------------------

def mergeListUnless(codingBoundsList, resolvedConflictsList, unUsableGeneNames, namesInResolvedConflictsList):
    
    for rec in codingBoundsList:
        
        #  If rec's name is in unUsable list or resolved list already move to next rec
        if rec[0] in unUsableGeneNames:
            continue
        elif rec[0] in namesInResolvedConflictsList:
            continue

        
        if boundaryRegion == 'up':
            if rec[3] in ['1','+']:
                resolvedConflictsList.append([rec[0], rec[1], rec[3], rec[8], rec[9], str(int(rec[9]) - int(rec[8]) + 1)])
            elif rec[3] in ['-1','-']:
                resolvedConflictsList.append([rec[0], rec[1], rec[3], rec[10], rec[11], str(int(rec[11]) - int(rec[10]) + 1)])
        
        elif boundaryRegion == 'down':
            if rec[3] in ['-1','-']:
                resolvedConflictsList.append([rec[0], rec[1], rec[3], rec[8], rec[9], str(int(rec[9]) - int(rec[8]) + 1)])
            elif rec[3] in ['1','+']:
                resolvedConflictsList.append([rec[0], rec[1], rec[3], rec[10], rec[11], str(int(rec[11]) - int(rec[10]) + 1)])
        else:
            print "WARNING: boundaryRegion variable should only be 'up' or 'down'.\nScript exiting."
            sys.exit()


#--------------------------------------------------       


# Strip trailing newlines
codingBoundsList = map(string.strip, codingBoundsList)
resolvedConflictsList = map(string.strip, resolvedConflictsList)
                            

# Convert these into lists of lists so that field vals can be interrogated and copied 
# Explode tab delimited strings of each record into list of values
JamesDefs.explodeDelimitedList(codingBoundsList, '\t')
JamesDefs.explodeDelimitedList(resolvedConflictsList, '\t')


len_codingBoundsList = len(codingBoundsList)
len_resolvedConflictsList = len(resolvedConflictsList)



# Populate unUsableList
unUsableGeneNames = []
i = 0
while i < len_resolvedConflictsList:
    
    if int(resolvedConflictsList[i][5]) < shortestUsableBdryReg:
        unUseableGene = resolvedConflictsList.pop(i)
        unUsableGeneNames.append(unUseableGene[0])
        unUseableGene = None
 
        #!!! when you pop _DO_NOT_ advance i and _DO_ reduce len_resolvedConflictsList to reflect new length
        len_resolvedConflictsList = len_resolvedConflictsList - 1 
    else:
        i=i+1
        

#  Create list of gene names still in resolvedConflictsList    
namesInResolvedConflictsList = []
for every in resolvedConflictsList:
    namesInResolvedConflictsList.append(every[0])
        
#  Append full boundary region records to resolved list unless name is in unUsableGeneNames
mergeListUnless(codingBoundsList, resolvedConflictsList, unUsableGeneNames, namesInResolvedConflictsList)


for each in resolvedConflictsList:
    reCombinedCoordsFILE.write('\t'.join(each)+'\n')


print 'Yay'

        
        
        
        
        
