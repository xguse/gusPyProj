#+++++++++++ Imports +++++++++++
import sys
import os
import optparse
import glob
from gusPyCode.defs.mpl_custom import autoPickColors
import random

import matplotlib as mpl
#mpl.use('Qt4Agg')
mpl.use('TkAgg')
from matplotlib import pylab as pl

from gusPyCode.defs.parseSCOPE import parseScopeXML
from gusPyCode.defs.JamesDefs import Bag


pl.rcParams['figure.figsize'] = 11, 8.5


#+++++++++++ Definitions +++++++++++

def random_color(color='#'):
    if len(color) == 7:
        return color
    else:
        color = color + hex(random.randrange(0, 255))[2:][:1]
        return random_color(color)

#colors = ['#FDFF00',
          #'#222222',
          #'#00FF00',
          #'#78006F',
          #'#00FFFF',
          #'#FF0000',
          #'#FFC300',
          #'#808080',
          #'#1B00FF',
          #'#FF00AE',]

def loadMotifs(args,sigLvl,filterMotifs):
    """load and store motifs from XML files giving each a unique ID and culling out the ones we want."""
    # --- load motifs ---
    motifDict = {}
    ID = 0 
    for path in args:
        tmpDict = parseScopeXML(path)
        # --- Give Each Motif its ID ---
        while tmpDict:
            k,v = tmpDict.popitem()
            # --- add the unique ID to the tuple key and the motif obj---
            # -- key --
            k = list(k) 
            k.insert(2,ID)
            k = tuple(k)
            # -- obj --
            v.motifID = ID
            # --- add motif to motifDict --
            motifDict[k] = v
            ID+=1
        
    # --- Cull Out The Ones We Want ---
    culledDict = {}
    while motifDict:
        k,v = motifDict.popitem()
        if float(v.sigvalue) < sigLvl:
            continue
        elif v.consensus in filterMotifs:
            culledDict[k] = v
    # --- Return Results ---
    for m in culledDict:
        print culledDict[m].consensus
    return culledDict
    
def groupMotifsByGene(geneList,motifDict):
    """group motifs by gene and sort by significance score. """ 
    motifsByGene = {}
    for g in geneList:
        genesMotifs = []
        # --- Collect Motifs in Gene ---
        for m in motifDict:
            if g in m[-1]:
                genesMotifs.append(motifDict[m])
        # --- Sort Motifs in Gene By SigVal ---
        genesMotifs.sort(key=lambda x: float(x.sigvalue))
        motifsByGene[g] = genesMotifs
    return motifsByGene
    

def loadGenes(geneListPath):
    """load and sort genes to be plotted alphabetically. """
    if not geneListPath:
        return None
    gList = map(lambda l: l.strip(), open(geneListPath, 'rU'))
    #gList.sort()
    gList.reverse()
    return gList

def plotGenes(ax,geneList,zorder):
    """create and add lines to represent genes and give them their zorder"""
    lines = []
    count = 0 
    for g in geneList:
        line = mpl.lines.Line2D((-2000,0),(count,count),color='black',linewidth=1, zorder=zorder)
        ax.add_line(line)
        lines.append(line)
        count+=1
        zorder+=1
    
    return lines


def plotMotifs4Gene(axes,patchDict,gene,motifsInGene,yVal,zorder):
    """create and add patches for a gene giving them their zorder """
    patches = []
    for motif in motifsInGene:
        for instance in motif.instances:
            if not instance[-1] == gene:
                continue
            color  = motif.color
            start  = min(int(instance[2]),int(instance[3]))
            width  = max(int(instance[2]),int(instance[3])) - min(int(instance[2]),int(instance[3])) + 1
            height = 0.5
            lowBnd = yVal-(float(height)/2)
            
            rect = mpl.patches.Rectangle( (start,lowBnd), color=color,
                                          width=width, height=height, zorder=None)
            
            patchDict[motif.consensus] = [rect] # stores last patch by consensus for legend
            patches.append(rect)
            zorder+=1
    for p in patches:
        ax.add_patch(p)
    
    
#def setMotifColors(motifDict,colors):
    #"""Assign colors to the filtered motifs in motifDict.  Adds
    #"self.color" attrib to the motif object directly."""
    #for m in motifDict:
        #motifDict[m].color = colors.pop()
        
def setMotifColorsManual(motifDict,colorDict):
    """Assign colors to the filtered motifs in motifDict.  Adds
    "self.color" attrib to the motif object directly."""
    print '-----'
    for motif in reversed(sorted(motifDict.values(),key=lambda x: float(x.sigvalue))):
        print '%s:%s' % (motif.consensus,colorDict[motif.consensus])
    for m in motifDict:
        motifDict[m].color = colorDict[motifDict[m].consensus]
        
def setMotifColorsRand(motifDict):
    """Build and use randomized colorDict to call
    setMotifColorsManual(motifDict,colorDict)."""
    colorDict = {}
    for m in motifDict.values():
        colorDict[m.consensus] = random_color()
    setMotifColorsManual(motifDict,colorDict)
    
def setMotifColors(motifDict):
    """AutoBuild and use colorDict to call 
    setMotifColorsManual(motifDict,colorDict)
    with contrasting colors."""
    if len(motifDict) < 4:
        colors = ['#0000FF','#008000','#FF0000']
    else:
        colors = list(autoPickColors(len(motifDict)))
    colorDict = {}
    for m in motifDict.values():
        colorDict[m.consensus] = colors.pop()
    setMotifColorsManual(motifDict,colorDict)
    
        
def buildLegend(axes,patchDict,consensusList):
    """Use motifDict and each motifs color to build a legend with consensus and color patch."""
    
    patches = []
    for c in consensusList:
        patches.append(patchDict[c])
    
    legend  = pl.legend(patches,consensusList, loc='upper left', numpoints=None,
                        markerscale=None, scatterpoints=3, scatteryoffsets=None,
                        prop=None, pad=None, labelsep=None, handlelen=None,
                        handletextsep=None, axespad=None, borderpad=None,
                        labelspacing=None, handlelength=None, handletextpad=None,
                        borderaxespad=None, columnspacing=None, ncol=1, mode=None,
                        fancybox=None, shadow=None, title=None,bbox_to_anchor=(1.01,0.99),
                        bbox_transform=None)

if __name__ == "__main__":
    #+++++++++++ Hard Code Figure Settings For Now +++++++++++
    title = ''
        
    print '\n\n\n'
    
    #+++++++++++ File Parseing Etc +++++++++++    
    usage = """python %prog [options]  scopeOut.xml [scopeOut.xml ...]"""
    parser = optparse.OptionParser(usage)
    
    parser.add_option('-g',dest="genes",type="str", default=None, 
                      help="""File listing genes to be ploted.  If None: plot all. (default=%default)""")
    parser.add_option('-o', dest="file_out", type='string',default='motifMap.pdf',
                      help="""File name for output (default=%default)""")
    parser.add_option('--sig', dest="sig_lvl", type='float',default=5.0,
                      help="""Only plot motifs with sig values >= this value. (default=%default)""")
    parser.add_option('--include', dest="include", type='str',default=None,
                      help="""REQUIRED: Comma separated list of consenus strings for motifs to plot. OPTIONALLY: Comma sep'd list of consensus:hexColorCode to manually set each motif color. (default=%default)""")
    parser.add_option('--rand-colors', dest="rand_colors", action="store_true",default=False,
                      help="""Choose colors for motifs randomly. (default=%default)""")
    defParams = '0.12,0.10,0.90,0.90,0.20,0.20'
    parser.add_option('--params', dest="params", type='str',default=defParams,
                      help="""Comma separated list of custom subplot parameters: 'left,bottom,right,top,wspace,hspace' (default=%default)""")
    
    
    (opts, args) = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help()
        exit(0)
    
    if len(args) == 0:
        print '\n*** ERROR: You must supply a at least one scopeXML file. ***\n'
        parser.print_help() 
        exit(1)
    
    # dumbass work-around for WING's sake (cant handle cmdLine glob patterns)
    if "*" in args[0]:
        args = glob.glob(args[0])   
    
    if opts.include:
        origInclude = opts.include
        if ":" in opts.include:
            tmpInclude = opts.include.split(',')
            opts.include = []
            colorDict = {}
            for elem in tmpInclude:
                elem = elem.split(':')
                colorDict[elem[0]] = elem[1]
                opts.include.append(elem[0])
        else:
            opts.include = opts.include.split(',')
    else:
        parser.print_help()
        print "\n*** ERROR: You must supply at least one value for '--include'. ***\n"
    
    if type(opts.params) == type(''):
        tmpParams = []
        for p in opts.params.split(','):
            tmpParams.append(float(p))
        if len(tmpParams) != 6:
            raise Exception('Malformed --params option. (%s)' % (opts.params))
        opts.params = tmpParams
        
    #+++++++++++ Load Motifs and GeneList +++++++++++
    motifDict    = loadMotifs(args,opts.sig_lvl,opts.include)
    geneList     = loadGenes(opts.genes)
    motifsByGene = groupMotifsByGene(geneList,motifDict)
    
    #+++++++++++ Init Figure +++++++++++
    fig = pl.figure()
    ax = fig.add_subplot(111)
    
    zorder = 1
    
    #+++++++++++ Plot Gene Lines +++++++++++
    geneLines = plotGenes(ax,geneList,zorder)
    pl.yticks(range(len(geneList)),geneList)
    
    #+++++++++++ Plot Motifs On Lines +++++++++++
    try:
        setMotifColorsManual(motifDict,colorDict)
    except NameError:
        if opts.rand_colors:
            setMotifColorsRand(motifDict)
        else:
            setMotifColors(motifDict)

        
    mPatches = []
    patchDict = {}
    for i in range(len(geneList)):
        mPatches.append(plotMotifs4Gene(ax,patchDict,geneList[i],motifsByGene[geneList[i]],i,zorder))
    
    #+++++++++++ Build Legend +++++++++++
    buildLegend(ax,patchDict,opts.include)
    
    None
    

    ax.autoscale_view()
    fig.subplots_adjust(left=opts.params[0],
                        bottom=opts.params[1],
                        right=opts.params[2],
                        top=opts.params[3],
                        wspace=opts.params[4],
                        hspace=opts.params[5])
    pl.savefig(opts.file_out)
    pl.show()

