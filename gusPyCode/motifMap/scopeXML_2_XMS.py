#!/usr/bin/env python
import optparse
import sys
from gusPyCode.defs.parseSCOPE import parseScopeXML,toXMSfile


if __name__ == "__main__":
    #+++++++++++ Hard Code Figure Settings For Now +++++++++++

        
    print '\n\n\n'
    
    #+++++++++++ File Parseing Etc +++++++++++    
    usage = """python %prog scopeOut.xml [scopeOut.xml ...]"""
    parser = optparse.OptionParser(usage)
        
    
    (opts, args) = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help()
        exit(0)
    
    for xml in args:
        out = xml.replace('.xml','.xms')
        motifs = parseScopeXML(xml)
        toXMSfile(motifs,out)
        