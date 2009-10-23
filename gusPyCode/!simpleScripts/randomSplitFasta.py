from TAMO.seq import Fasta


fasFile = '/Users/biggus/Documents/James/Writings_Talks/Grants/09_Feb/PrelimData_Grant_Feb09/Clus2_247genes.fas'
oFile1= '/Users/biggus/Documents/James/Writings_Talks/Grants/09_Feb/PrelimData_Grant_Feb09/Clus2_247genes.sample2.fas'
oFile2= '/Users/biggus/Documents/James/Writings_Talks/Grants/09_Feb/PrelimData_Grant_Feb09/Clus2_247genes.test2.fas'

firstDic, secDic = Fasta.random_split(fasFile,0.25)

Fasta.write(firstDic,oFile1)
Fasta.write(secDic,oFile2)

print 'done'