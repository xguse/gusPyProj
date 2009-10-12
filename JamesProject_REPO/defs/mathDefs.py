import numpy as np
import scipy as sp

def mean(listOfNums):
    tot = 0
    for each in listOfNums: 
        tot += each
    return float(tot)/len(listOfNums)

def line(m,b):
    """Given m(slope) and b(y-intercept) return y"""
    

def medianAbsDev(dataVector):
    """
    Takes a 1-D data vector.
    Returns (MAD, median).
    """
    median  = np.median(dataVector)
    absDevs = [np.abs(x-median) for x in dataVector]
    MAD = np.median(absDevs)
    
    return (MAD, median)

def stdDv(dataVector,kind='mean',df=0): 
    """
    stdDv(dataVector,type='mean',df=0)
    df = how many degOfFreedom are to subtracted from sample number
    Returns (StdDv, medianOrMean)
    """
    assert kind in ['median','mean'], \
           "ERROR[stdDv(dataVector,kind='mean',df=0)]: Valid values for 'kind' are 'mean' or 'median'"
    
    if kind == 'mean':
        m  = np.mean(dataVector)
    else:
        m  = np.median(dataVector)
        
    variance  = np.sum(np.square([x-m for x in dataVector]))/(len(dataVector)-df)
    stDev     = np.sqrt(variance)
    
    return (stDev,m)
    






changeLog = """2009-08-29 -- module created.
2009-08-29 -- added medianAbsDev(dataVector)
2009-08-29 -- added stdDv(dataVector,kind='mean',df=0)
2009-09-14 -- added mean(listOfNums)
"""
