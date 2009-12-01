from math import pow

def qual2prob(qualChar):
    """Takes a single Phred quality score character and returns a the probability
    that the corresponding base was called incorrectly.
    Q_phred = -10 * log10(p)
    p = 10^(Q_phred/-10)
    """
    q = ord(qualChar)-64
    return pow(10,(float(q)/-10))