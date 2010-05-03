import optparse
import sys
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


(opts, args) = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    exit(0)
if (not opts.fasta) or (not opts.regs) or (not opts.out):
    parser.print_help()
    print '\n'
    raise Exception, 'ERROR: Missing a required input.'






optionsObj = {'kmerSize':opts.ksize,'kmerRange':opts.krange} # need to change the use of this to be how mortals think of numbers
posArgs    = [opts.fasta,opts.regs]
outFile    = opts.out

tWrap = TamoWrap(optionsObj,posArgs)
tWrap.go()

#for l in tWrap.toFile:
    #print l.strip()

print 'num of keepers = %s' % (len(tWrap.output))

outFile = open(outFile, 'w')
outFile.writelines(tWrap.toFile)