import sys

def statusBar(numOfActionsAccomplished,numOfActionsTotal):
    """
    Takes numberOfActionsAccomplished and numberOfActionsTotal and writes a staus bar to indicate
    the progress.
    """
    
    blocks = int(numOfActionsAccomplished * 40 / numOfActionsTotal)
    #print 'blocks: %s    40-blocks: %s' % (blocks,40-blocks)
    
    sys.stdout.write('[%s%s]\r' % ('='*blocks, ' '*(39-blocks)))
    sys.stdout.flush()
    sys.stdout.write("%s" % ('\b'*42))
    sys.stdout.flush()

tot = 86410

for i in range(1,tot):
    blocks = int(i * 40 / tot)
    #print 'blocks: %s    40-blocks: %s' % (blocks,40-blocks)

    sys.stdout.write('[%s%s]' % ('='*blocks, ' '*(39-blocks)))
    sys.stdout.flush()
    sys.stdout.write("%s" % ('\b'*42))



