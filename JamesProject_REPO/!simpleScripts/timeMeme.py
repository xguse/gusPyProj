from time import time
import os

cmd = "meme Clus2_247genes.fas -dna -revcomp -minw 6 -maxw 12 -mod zoops  -nostatus -text -maxsize 426302 -p > testMemeTime.meme.txt"

t1 = time()
FID = os.popen(cmd ,'r')
t2 = time()

print