from gusPyCode.defs import JamesDefs



coordList =  open('/Users/biggus/Documents/MBGB/Rotations/James/Data/2KB/2kb_Sequence/2kb_Aedes/2kb_aedesProcessed/2KBup_aedesUsableBdryCoordsTSS.txt','rU').readlines()
listList =  open('/Users/biggus/Desktop/testAPI/wtfAPI/AedesMissingEnsemblList.txt','rU').readlines()

outFile = open('/Users/biggus/Desktop/testAPI/wtfAPI/AedesMissingEnsemblListCoords.txt','w')

rList = []

for line in coordList:
    fields = line.split('\t')
    check = fields[0]+'\n'
    if check in listList:
        rList.append(line)

outFile.writelines(rList)