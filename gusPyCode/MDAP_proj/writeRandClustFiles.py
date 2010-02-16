from gusPyCode.MDAP_proj.MDAP_defs import genRandClusters
from gusPyCode.defs.JamesDefs import mkdirp
from TAMO.seq import Fasta
from numpy import average
import pprint
import glob
import optparse



# ----- Parse CmdLine -----
print '\n\n\n'

    
usage = """python %prog [options]  geneListFile fastaFile"""
parser = optparse.OptionParser(usage)

parser.add_option('-N',dest="N",type="int", default=1, 
                  help="""Number of ctrl sets to make (default=%default)""")
parser.add_option('-o', dest="out_dir", type='string',default='./geneNameCtrlClusters',
                  help="""Dir for output (default=%default)""")
parser.add_option('-f', dest="make_fasta", action="store_true",default=False,
                  help="""Produce relavent fasta files too. (default=%default)""")


(opts, args) = parser.parse_args()

# --- A Little Extra Input Validation ---
if len(args) < 2:
    parser.print_help()
    print '\nERROR: Both geneListFile and fastaFile are required!'
    exit(1)



geneNames = map(lambda l: l.strip(),open(args[0], 'rU').readlines())
totalSeqs = Fasta.file2dict(args[1])
randClusterLists = genRandClusters(geneNames,totalSeqs,N=opts.N, keepLen=1)

# -- Make Out Folder --
mkdirp(opts.out_dir)
    

for i in range(len(randClusterLists)):
    oFileName = args[0].replace('.txt','randomGeneNames_%s.txt' % (i)).split('/')[-1]
    oFile = open('%s/%s' % (opts.out_dir,oFileName), 'w')
    for name in randClusterLists[i]:
        oFile.write(name+'\n')
    oFile.close()
    # --- If Asked, Create Fastas ---
    if opts.make_fasta:
        fNames  = map(lambda l: l.strip(),open('%s/%s' % (opts.out_dir,oFileName), 'rU').readlines())
        fastas = {}
        for f in fNames:
            fastas[f] = totalSeqs[f]
        Fasta.write(fastas,'%s/%s' % (opts.out_dir,oFileName.replace('.txt','.fas')))
    del(oFile)
    
print 'Done.'


            
print "Original list:"
origLens = []
for g in geneNames:
    origLens.append(len(totalSeqs[g]))
origLens.sort()
print origLens
print 'Avg: %.3f' % (average(origLens))
    
print "Randomized lists:"
for r in randClusterLists:
    randLens = []
    for g in r:
        randLens.append(len(totalSeqs[g]))
    randLens.sort()
    print randLens
    print 'Avg: %.3f' % (average(randLens))
    