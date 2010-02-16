import optparse

print '\n\n'


# --------- Process User Input --------- 
usage = """python %prog [options]  inFile [xtraFile xtraFile]"""
parser = optparse.OptionParser(usage)

parser.add_option('-o', dest="out_file", type="string", default=True,
                  help="""FileName for outFile""")

(opts, args) = parser.parse_args()


if len(args) < 2:
    parser.print_help()
    print '\nERROR: Please supply at least 2 inFiles.'
    exit(1)

if opts.out_file == True:
    parser.print_help()
    print '\nERROR: Please supply a fileName for the outFile: "-o fileName4out".'
    exit(1)

for f in args:
    if opts.out_file.split('/')[-1] == f.split('/')[-1]:
        print "ERROR: You supplied the name of one of your inFiles for the outFile.  This is not allowed."
        exit(1)
    
    
# --------- Append Files ---------

fOut = open(opts.out_file,'w')

for inFile in args:
    print "Appending %s..." % (inFile)
    for line in open(inFile,'rU'):
        fOut.write(line)
        
print "Done"
