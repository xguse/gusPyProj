"""Convert a Transfac file into a possum matrix.

Allows saving to separate files.

Implemented by Kenny Daily 2009/2010
Institute of Genomics and Bioinformatics
University of California, Irvine

All software and related material contained 
herein can be downloaded freely for academic,
non-commercial, research use only.

For any other use, please contact Pierre Baldi at pfbaldi _at_ ics uci edu. 

Copyright 2009-2010

Contact: kdaily _at_ ics uci edu

"""

import sys
import os
import optparse

from gusPyCode.defs import motifparsers

def main():
    """convert a Transfac file into a possum matrix/matrices.

    Allows saving to separate files.
    
    """

    usage = "%prog [options] transfac_file [> possum_file]" 
    parser = optparse.OptionParser(usage)

    group = optparse.OptionGroup(parser, "File split options", "Options for splitting the matrices outputted to separate files")
    group.add_option("-s", "--split", action="store_true", dest="split", default=False,
                      help="Split the file by accession number. If given, will write the files using the name in the AC field, appended with '.possummatrix',  to the directory given by the -o flag. Default is to output a single file.")
    group.add_option("-o", "--output_directory", type="string", dest="output_dir", default='.',
                      help="Used only with the -s (--split) flag to set the output directory.")
    parser.add_option_group(group)

    (opts, args) = parser.parse_args()
    
    if len(args) == 0:
        parser.print_help() 
        exit(1)
        
    transfac_file = args[0]

    ## load pwm
    transfac_data = motifparsers.dbparser['transfac'](transfac_file)
        
    if opts.split:
        # save out a separate file for each accession number
        for acc in transfac_data.getMotifList():
            outfile = file(os.path.join(opts.output_dir, '%s.possummatrix' % acc), 'w')
            outfile.write(transfac_data.exportDataAsPossum(acc))
            outfile.close()
    else:
        outfile = sys.stdout
        outfile.write(transfac_data.exportDataAsPossum())
        outfile.close()
    
if __name__ == '__main__':
    main()
