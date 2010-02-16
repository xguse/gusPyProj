#!/usr/bin/env python
"""Build an annotationDB from the given gff file

"""

import optparse
import sys

from pygr import worldbase
from pygr import annotation
from pygr import cnestedlist
from pygr import metabase
import sqlite3

from gusPyCode.pygr.gffparser import read_for_pygr
from gusPyCode.pygr.gPygr import simpleGFF2PygrSQLite,getTableNames
from BCBio.GFF import GFFParser


def main():
    """Build an annotation from the given gff file
    """
    
    usage = """Build and save the annotations defined in the given gff files
    Saves an annotationDB (representing the file itself) and creates a mapping 
    in the form genome[chromosome][100:200].officialGenes"""
    parser = optparse.OptionParser("%prog [options] data1.gff [data2.gff ...]\n"+usage)
    parser.add_option("--genome_resource", '-g', dest="genome_resource", type="string",
                      help="""The pygr resource for the genome, eg, 'Bio.Seq.Genome.TRICA.triCas3'""")
    #parser.add_option("--annotationDB_resource", '-a', dest="annotationDB_resource", type="string",
                      #help="""Where to save the created annotationDB. eg, 
                      #Bio.Annotation.TRICA.triCas3.officialGenes""")
    parser.add_option("--sqlDB_resource", '-s', dest="sqlDB_resource", type="string",
                      help="""Where to save the created sqlDB and a unique file name eg, 
                      Bio.Annotation.TRICA.triCas3.features_sqlDB,gffDB_v1""")
    parser.add_option("--save_pathstem", '-p', dest="pathstem", type="string", 
                      help="""The file to save the resource to, eg,
                    '/home/baldig/projects/genomics/pygrdata/annotations/fly/triCas3_official_genes'""")
    parser.add_option("--map_resource", '-m', dest="map_resource", type="string",
                      help="""the resource to save the annotationDB->Genome map,
                      saved both to worldbase and to worldbase.schema, eg,
                      'Bio.Annotation.TRICA.triCas3.BeetleBase.officialGenesMap""")
    parser.add_option("--bind_attribute", '-b', dest="bind_attribute", type="string", 
                      help="""The attribute to access annotationDB from genome region, eg, 
                      'officialGenes' would be accessible via triCas3['ChLG2'][100:200].officialGenes 
                      Default is not to bind an attribute to genome""")


    (opts, args) = parser.parse_args()

    if len(args) < 1: 
        parser.print_help()
        print 'Please specify at least one gff file to read'
        sys.exit(-1)
    if None in [opts.genome_resource, opts.pathstem, opts.map_resource]:
        parser.print_help()
        print 'Required options: genome_resource, sqlDB_resource, pathstem, map_resource'
        sys.exit(-1)
    if opts.sqlDB_resource.count(',') != 1:
        parser.print_help()
        print 'Error: sqlDB_resource must be comma separated string with exactly one comma.'
    else:
        opts.sqlDB_resource = opts.sqlDB_resource.split(',')
    try :
        w = worldbase(opts.sqlDB_resource[0])
        parser.print_help()
        print "Warning: sqlDB_resource already exists.  Please select a new name."
        exit(-1)
    except WorldbaseNotFoundError:
        pass
    
    
    print '# Loading original genome db'
    genome = worldbase(opts.genome_resource)
    #annotDB = annotation.AnnotationDB(None, genome, opts.bind_attribute, 
                                        #filename=opts.pathstem + '_annotDB', mode='c', verbose=False)
    sqlDB    = sqlgraph.SQLiteServerInfo('%s/%s.sqlite' %(opts.pathstem,opts.sqlDB_resource[1]))
    gff2lite = simpleGFF2PygrSQLite(sqlDB)
    nlmsa    = cnestedlist.NLMSA(opts.pathstem, 'w', pairwiseMode=True, bidirectional=False)
    
    
    for filename in args:
        print '# adding to sqlDB from %s' % filename
        gff2lite.update(filename)
    
    tableNames = gff2lite.getTableNames()
    for table in tableNames:
        
    
    
        
    #for row in read_for_pygr(fileIn):
        #curAnnot = annotDB.new_annotation(index, row)
        #nlmsa.addAnnotation(curAnnot)
        #index += 1
    #annotDB.close() # Flush annotation data to disk
    
    print '# building NLMSA from all gff files'
    nlmsa.build(saveSeqDict=True)
    print '# saving annotationDB and NLMSA to worldbase as %s and %s' % (opts.annotationDB_resource,
                                                                        opts.map_resource)
    annotDB.__doc__ = 'Combined gff annotationDB from files %s on genome %s' % (', '.join(args), 
                                                                                opts.genome_resource)
    nlmsa.__doc__ = 'Mapping of %s, from gff files %s onto genome %s' % (opts.annotationDB_resource,
                                                                            ', '.join(args),
                                                                            opts.genome_resource)
    worldbase.add_resource(opts.annotationDB_resource, annotDB)
    worldbase.add_resource(opts.map_resource, nlmsa)

    if opts.bind_attribute:
        print '# saving worldbase schema with bindAttrs=(%s)' % opts.bind_attribute
        genome_annotDB_relation = metabase.ManyToManyRelation(genome, annotDB, bindAttrs=(opts.bind_attribute,))
        genome_annotDB_relation.__doc__ = 'GFF based mapping from %s to genome %s' % (opts.annotationDB_resource,
                                                                                        opts.genome_resource)
        worldbase.add_schema('%s' % opts.map_resource, genome_annotDB_relation)
                                
    
    print '# committing worldbase resources'
    worldbase.commit()

if __name__ == "__main__":
    main()
    
