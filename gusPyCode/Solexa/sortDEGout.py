import optparse

print '\n\n'


# --------- Process User Input --------- 
usage = """python %prog inFile [xtraFile xtraFile]"""
parser = optparse.OptionParser(usage)



(opts, args) = parser.parse_args()


if len(args) < 1:
    parser.print_help()
    print '\nERROR: Please supply at least 1 inFiles.'
    exit(1)


# --------- Iterate Through Files ---------
for inFile in args:
    data = []
    print "Sorting %s..." % (inFile)
    # --------- Collect Lines ---------
    for line in open(inFile,'rU'):
        if "GeneNames" in line:
            continue
        l = line.split('\t')
        if l[3] != 'NA':
            data.append(l)
    # --------- Sort Lines ---------
    data.sort(key=lambda x: abs(float(x[3])))
    data.reverse()
    # --------- Write Lines ---------
    fOut = open(inFile+".sorted.txt",'w')
    for item in data:
        fOut.write("%s" % ('\t'.join(item)))
    fOut.close()
        

        
print "Done"
