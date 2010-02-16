import matplotlib as mpl
from matplotlib import pylab as pl


class Seq(object):
    """Class to manage motifs found in a gene.
    """
    def __init__(self,motifData):
        '''motifData = list(motifID, strand<-/+>,start,end)
        '''
        self.motifs = motifData