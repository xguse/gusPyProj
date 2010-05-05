import optparse
import sys
from warnings import warn
from gusPyCode.MDAP_proj.MD_wrappers import TamoWrap

# TamoWrap -> def __init__(self, optionsObj, posArgs):

#+++++++++++ File Parseing Etc +++++++++++    
usage = """python %prog [options] """
parser = optparse.OptionParser(usage)

parser.add_option('--kSize',dest="ksize",type="int", default=7, 
                  help="""k-mer size. (default=%default)""")
parser.add_option('--kRange', dest="krange", type='string',default='6,10',
                  help="""k-mer range (python range-type). (default=%default)""")
parser.add_option('-o', dest="out", type='string',default=None,
                  help="""REQUIRED: File name for output (default=%default)""")
parser.add_option('--fas', dest="fasta", type='str',default=None,
                  help="""REQUIRED: Fasta file of ALL regulatory sequences. (default=%default)""")
parser.add_option('--regs', dest="regs", type='str',default=None,
                  help="""REQUIRED: File with list of regulated genes. (default=%default)""")
parser.add_option('--trim', dest='trim', type='str',default=None, \
				  help="""Define substring of all sequences to be used.
                  Exp: "-2000,None" -> last 2000 bp of all sequences. (default=%default)""")


(opts, args) = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    exit(0)
if (not opts.fasta) or (not opts.regs) or (not opts.out):
    parser.print_help()
    print '\n'
    raise Exception, 'ERROR: Missing a required input.'
if opts.trim:
    opts.trim = opts.trim.split(',')
    if opts.trim[0] == 'None':
        opts.trim[0] = None
    else:
        opts.trim[0] = int(opts.trim[0])
        
    if opts.trim[1] == 'None':
        opts.trim[1] = None
    else:
        opts.trim[1] = int(opts.trim[1])





optionsObj = {'kmerSize':opts.ksize,'kmerRange':opts.krange} # need to change the use of this to be how mortals think of numbers
posArgs    = [opts.fasta,opts.regs]
outFile    = opts.out

tWrap = TamoWrap(optionsObj,posArgs)

if opts.trim:
    shorties = 0
    for i in tWrap.allSeqs.probes:
        oLen = len(tWrap.allSeqs.probes[i])
        tWrap.allSeqs.probes[i] = tWrap.allSeqs.probes[i][opts.trim[0]:opts.trim[1]]
        nLen = len(tWrap.allSeqs.probes[i])
        #print 'old:%s, new:%s' % (oLen,nLen)
        if not nLen < oLen:
            shorties+=1
if shorties/float(len(tWrap.allSeqs.probes)) > 0.25:
    warn("""WARNING: more than 1/4 of the total sequences were shorter or equal to the length of the requested substring.""")
tWrap.linkedSeqs_seqs = tWrap.allSeqs.seqs_from_ids(tWrap.linkedSeqs_ids)
tWrap.go()

#for l in tWrap.toFile:
    #print l.strip()

print 'num of keepers = %s' % (len(tWrap.output))

outFile = open(outFile, 'w')
outFile.writelines(tWrap.toFile)