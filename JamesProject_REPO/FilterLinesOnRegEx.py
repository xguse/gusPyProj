#!/usr/bin/env python

import sys
import re


if len(sys.argv) != 3:
    print "Usage: python tab-ID_filter.py comparison_file filter_ID_file"
    sys.exit(1)

comparison_file = sys.argv[1]
filter_file =  sys.argv[2]

comparison_List = map(lambda line : line.strip(), open(comparison_file, 'rU').readlines())
filter_List     = map(lambda line : line.strip(), open(filter_file, 'rU').readlines())



for m in filter_List:
    # add commas to both sides of mottif for boundings
    m_Bounds = ','+m+','
    m_RegEx = re.compile(m_Bounds, re.IGNORECASE)
    
    hitList = []
    for line in comparison_List:
        if m_RegEx.search(line):
            hitList.append(line)
    print "lines matching %s exactly:" % (m)
    for each in hitList:
        print each
        
    print '-------------------------------'
    
    