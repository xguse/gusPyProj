

#---------

def simulate_sequence(length):
    import random

    dna = ['A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 'S', 'W']  
    sequence = ''  
    for i in range(length):  
        sequence += random.choice(dna)

    totDegen = countDegenInSeq(sequence)

    if totDegen > 2:
        sequence = None
    
    return sequence

def countDegenInSeq(seq):
    degenSet = ['R', 'Y', 'M', 'K', 'S', 'W']

    explodedSeq = list(seq)

    totalDegen = 0

    for each in explodedSeq:
        if each in degenSet:
            totalDegen+=1

    return totalDegen

#---------


#========================= User Defined Variables =========================
InFile  = ''
KmerLength = 8
numOfRand_Kmers = 75

outFile = ''

#==========================================================================

#InFile = map(lambda line: line.strip(), open(mofifInFile, 'rU').readlines())


rand_Kmers = []


while len(rand_Kmers) <= numOfRand_Kmers:
    newRand_Kmer = simulate_sequence(KmerLength)
    if newRand_Kmer:
        rand_Kmers.append(newRand_Kmer)
        
for each in rand_Kmers:
    print each
    