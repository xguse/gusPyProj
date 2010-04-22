aCnts = {'1':{3:4,
              4:3,
              5:5,
              6:7,
              10:10,
              11:9,}}

def _writeWiggle(trackName, trackDescription, allCounts, wigOut, win=1):
    """Writes the allCounts dictionary out to a wiggle file."""
    #wigFile = open(wigOut, "w")
    #wigFile.write("track type=wiggle_0 name='%s' description='%s' visibility=2\n" % (trackName,
                                                                                     #trackDescription))
    print "track type=wiggle_0 name='%s' description='%s' visibility=2\n" % (trackName,
                                                                             trackDescription)
    assert type(win) == type(1), \
           'Error: win must be integer. You gave: %s' % (win)

    for name in allCounts.keys():
        start = 0
        end = max(allCounts[name].keys())

        #wigFile.write("fixedStep chrom=%s start=%s step=%s span=%s\n" % (name, start+1, win, win))
        print "fixedStep chrom=%s start=%s step=%s span=%s\n" % (name, start+1, win, win)
        i=0
        #for i in range(start, end):
        while i <= end:
            winCounts = []
            for pos in range(i,i+win):
                if allCounts[name].has_key(pos):
                    winCounts.append(allCounts[name][pos])
                else:
                    winCounts.append(0)
            winAvg = sum(winCounts)/float(len(winCounts))
            #wigFile.write("%.4f\n" % (winAvg))
            print "%.4f\n" % (winAvg)
            i+=win

    #wigFile.close()

    
def _writeWiggleVar(trackName, trackDescription, allCounts, wigOut):
    """Writes the allCounts dictionary out to a wiggle file."""
    wigFile = open(wigOut, "w")
    wigFile.write("track type=wiggle_0 name='%s' description='%s' visibility=2\n" % (trackName,
                                                                                     trackDescription))

    for name in allCounts.keys():
        start = 0
        end = max(allCounts[name].keys())

        wigFile.write("variableStep chrom=%s span=1\n" % (name))

        for pos in sorted(allCounts[name].keys()):
            wigFile.write("%s\t%s\n" % (pos,allCounts[name][pos]))


    wigFile.close()