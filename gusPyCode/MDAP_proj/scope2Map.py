import sys
import os
import optparse
import glob

import weblogolib as wll
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram as gd
from reportlab.lib.units import cm
from gusPyCode.defs.parseSCOPE import parseScopeXML
#import cPickle
from gusPyCode.defs.JamesDefs import Bag


def mkdirp(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)
        
def testForGenes(geneList, motif):
    geneList = set(geneList)
    hitGenes = set(motif[2])
    if len(geneList) == len(geneList.intersection(hitGenes)):
        return True
    else:
        return False
    
def makeLogo(motifList,outFile,oType='pdf'):
    assert oType in ['pdf','eps'],\
           'ERROR: oType value incorrect. Please one of: %s' % (['pdf','eps'])
    seqs = motifList
    for i in range(len(seqs)):
        seqs[i] = wll.Seq(seqs[i])
        
    seqs = wll.SeqList(seqs)
    seqs.alphabet = wll.which_alphabet(seqs)
    
    data = wll.LogoData.from_seqs(seqs)
    options = wll.LogoOptions() #color_scheme='classic'
    #options.title = 'A Logo Title'
    format = wll.LogoFormat(data, options)
    fout = open('%s.%s' % (outFile,oType), 'w')
    if oType == 'eps':
        wll.eps_formatter( data, format, fout)
    else:
        wll.pdf_formatter( data, format, fout)

def mInst2SeqFeature(mtfInst,prmtrLen):
    """Take a motif instance and return a SeqFeature with location obj.
    """

    if mtfInst.strand == '-':
        strand = -1
        start  = prmtrLen+int(mtfInst.end)
        stop   = prmtrLen+int(mtfInst.begin)+1
    elif mtfInst.strand == '+':
        strand= 1
        start  = prmtrLen+int(mtfInst.begin)
        stop   = prmtrLen+int(mtfInst.end)+1
    else:
        print 'Motif instance.strand was not "-" or "+".'
        exit(1)
    seqFeat = SeqFeature(FeatureLocation(start,stop),strand=strand)
    return seqFeat

    

#usage = '\n\npython %s outDir "geneName geneName ...;query" scopeOut.xml [scopeOut.xml]' % (sys.argv[0].split('/')[-1])

#if len(sys.argv) < 3:
    #print usage
    #exit(0)

# ----- Parse CmdLine -----
print '\n\n\n'

    
usage = """python %prog [options]  scopeOut.xml [scopeOut.xml ...]"""
parser = optparse.OptionParser(usage)

parser.add_option('-p',dest="promoterLen",type="int", default='2000', 
                  help="""Length of promoters: (default=%default)""")
parser.add_option('-o', dest="outDir", type='string',default='./scopeFigs',
                  help="""Dir for output (default=%default)""")
parser.add_option('-q',dest="geneSearch",type="string", default=None, 
                  help="""REQUIRED:Quoted string with required genes\nbefore a ';' and non-essential after.\n "reqGene reqGene;nReqGene" """)



(opts, args) = parser.parse_args()

# dumbass work-around for WING's sake (cant handle cmdLine glob patterns)
if "*" in args[0]:
    args = glob.glob(args[0])

if not opts.geneSearch:
    print '\n*** ERROR: You must supply a geneSearch string. ***\n'
    parser.print_help() 
    exit(1)

# Remove trailing '/' if present.
if opts.outDir[-1] == '/':
    opts.outDir = outDir[:-1]
# -------------------------
    
    

    
# 1. Open all files into own scopeOut dicts
inDict = {}
for path in args:
    inDict[path] = parseScopeXML(path)
    
# Parse GeneNames
if opts.geneSearch:
    if ';' in opts.geneSearch:
        geneNames,query = opts.geneSearch.split(';')
        geneNames = geneNames.split() 
    else:
        geneNames = opts.geneSearch.split()
        query = None
#else:
    

# Print list of motifs and File names that hit genes listed:
hits = {}
for f in inDict:
    #tmp = {}
    for motif in inDict[f]:
        if testForGenes(geneNames,motif):
            hits[(f,motif)] = inDict[f][motif]

hitKeys = hits.keys()
hitKeys.sort(key=lambda x: (x[0],x[1]))

# ------- Create Dir For Output -------
mkdirp(opts.outDir)


# ------- Format and Write Output -------
if hitKeys:
    for k in hitKeys:
        if float(hits[k].sigvalue) < 5:
            del hits[k] # Remove from hits if SigVal less than 5.
        else:
            pass
            #fName = k[0]
            #if '/' in fName:
                #fName = fName.split('/')[-1]
            #queries = []
            #for g in geneNames+[query]:
                #if g in hits[k].genes:
                    #queries.append(g)
            #queries.sort()
            #print '%s\t%s\t%s\t%s\t%s\t%s\t%s' % (fName,
                                              #hits[k].sigvalue,
                                              #hits[k].sigValRank,
                                              #hits[k].consensus,
                                              #round(len(geneNames)/float(len(hits[k].genes))*100,1),
                                              #hits[k].algorithm,
                                              #'_'.join(queries))
            #motifList = [x.sequence for x in hits[k].instances]
            #makeLogo(motifList,
                     #'%s/%s.%s_%s'% (opts.outDir,fName,hits[k].consensus,hits[k].sigValRank),
                     #oType='pdf')
            #makeLogo(motifList,
                     #'%s/%s.%s_%s'% (opts.outDir,fName,hits[k].consensus,hits[k].sigValRank),
                     #oType='eps')
    
else:
    print 'No motifs matched all of your genes. :('
    exit(0)
    
    


# ----- Cycle through genes of interest -----
for gene in geneNames+[query]:
    print 'Mapping motifs for %s...' % (gene)
    diag  = gd.Diagram() # (1)
    # ----- Build feature sets -----
    for k in hits:
        # decide if motif k hits our gene
        if gene in hits[k].genes:
            # If so, open new trk and fSet
            newTrack = diag.new_track(1,name=hits[k].consensus,greytrack=1,greytrack_labels=1,scale=1) # (2)
            newSet   = newTrack.new_set()
            for instance in hits[k].instances:
                # Is this instance of motif k in our current gene?
                if instance.gene == gene:
                    # if so, create feature and add it to fSet
                    sFtr = mInst2SeqFeature(instance,opts.promoterLen)
                    newSet.add_feature(sFtr,color=colors.red, label=True)
            
        # ----- Draw/Save Maps -----
        ##diag.draw(format='linear', pagesize=(2*cm,19*cm), orientation='landscape', fragments=1, start=0, end=2000)
        diag.draw(format='linear', orientation='portrait', fragments=1, start=0, end=opts.promoterLen)
        diag.write("%s/%s.pdf" % (opts.outDir,'_'.join([gene,'motifMaps'])), "pdf")
print 'Done.'
