import sys
import optparse
from xml.dom import minidom

def getContig(xmlDoc):
    return str(dasDoc.getElementsByTagName('SEGMENT')[0].getAttribute('id'))

def getFeatID(feature):
    return str(feature.getAttribute('id'))
def getFeatMethod(feature):
    return str(feature.getElementsByTagName('METHOD')[0].firstChild.data)
def getFeatStart(feature):
    return str(feature.getElementsByTagName('START')[0].firstChild.data)
def getFeatEnd(feature):
    return str(feature.getElementsByTagName('END')[0].firstChild.data)
def getFeatOri(feature):
    return str(feature.getElementsByTagName('ORIENTATION')[0].firstChild.data)
def getFeatGroup(feature):
    return str(feature.getElementsByTagName('GROUP')[0].getAttribute('id'))
def getFeatGroupType(feature):
    return str(feature.getElementsByTagName('GROUP')[0].getAttribute('type'))


if __name__ == "__main__":
    
    
    #+++++++++++ File Parseing Etc +++++++++++    
    usage = """python %prog inFile [options]"""
    parser = optparse.OptionParser(usage)
    parser.add_option('-o',dest="out_file",type="str", default="outfile.txt", 
                      help="""Space-less string for the display name for this track. (default=%default)""")

    
    (opts, args) = parser.parse_args()
    
    if len(args) < 1:
        parser.print_help()
        exit()
    
    
    outData = []
    for path in args:
        # +++ create DOM +++
        dasDoc = minidom.parse(path)
        contig = getContig(dasDoc) 
        for feature in dasDoc.getElementsByTagName('FEATURE'):
            feat_id = getFeatID(feature)
            method  = getFeatMethod(feature)
            start   = getFeatStart(feature)
            end     = getFeatEnd(feature)
            ori     = getFeatOri(feature)
            group   = getFeatGroup(feature)
            grpTyp  = getFeatGroupType(feature)
            
            outData.append([contig,
                           group,
                           feat_id,
                           start,
                           end,
                           ori,
                           method,
                           grpTyp])

    
    outFile = open(opts.out_file, 'w')
    outFile.write('contig\tgroup\tfeat_id\tstart\tend\tori\tmethod\tgrpTyp\n')
    for line in outData:
        outFile.write('%s\n' % ('\t'.join(line)))
    outFile.close()
            