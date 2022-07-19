import sys, os
from math import log

def H(n):
    """
    Found this answer on: https://stackoverflow.com/questions/404346/python-program-to-calculate-harmonic-series
    Returns an approximate value of n-th harmonic number.
    http://en.wikipedia.org/wiki/Harmonic_number
    """
    # Euler-Mascheroni constant
    gamma = 0.57721566490153286060651209008240243104215933593992
    return gamma + log(n) + 0.5/n - 1./(12*n**2) + 1./(120*n**4)

if __name__=="__main__":
    
    numSegregatingSites=int(sys.argv[1])
    numIndividuals=int(sys.argv[2])
    len_chr=int(sys.argv[3])
    numHaplotypes=numIndividuals*2
    a_n=H(numHaplotypes-1)
    theta=numSegregatingSites/(a_n*len_chr)
    print(theta)

