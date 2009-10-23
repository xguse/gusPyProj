import random
import sys
#========================= User Defined Variables =========================
list = sys.argv[1]

sampleSize = int(sys.argv[2])

#==========================================================================


list = map(lambda line: line.strip(), open(list, 'rU').readlines())

sampleList = []
for i in range(0,sampleSize):
    sampleList.append(random.choice(list))
    
    
for each in sampleList:
    print each