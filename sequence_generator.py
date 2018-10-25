#!/usr/bin/env python

"""
This module generates random numbers from a given function"s
"""
__author__ = "Pratha Sah and Shweta Bansal"
__copyright__ = "Copyright (C) 2013 Pratha Sah"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Pratha Sah"
__email__ = "ps875@georgetown.edu"

import numpy as np
import math
import random as rnd

###################################################################################################

def regular_sequence (N, mean, seqtype):
	"""returns a sequence with N entries and values symmetric around the mean"""
        seq=[int(round(mean))]*(N)
        return seq
        
##############################################################################################   
def poisson_sequence (N, mean, seqtype):
	"""returns a sequence with N entries and values following a Poisson distribution with mean = mean"""
	seq= list(np.random.poisson(lam=mean, size=(N)))
        return seq

##############################################################################################   

def scalefree_sequence (N, mean, seqtype):
	"""returns a sequence with N entries and values following a Power-law distribution with mean = mean"""
	if seqtype =="modulesize": alpha=10
	else: alpha=1/(1.0*mean)*10
        condition=False
        tol= 5 # set initial tolerance high to enter the loop
        while tol> 0.2:
            #print alpha,
            seq1= list((np.random. power(alpha, size=(N))))
            seq1=[abs(1-x) for x in seq1]
            if seqtype == "modulesize": seq=[int(num*(N*mean-1)) for num in seq1]
            if seqtype == "simple_degree": seq=[int(num*(N-2))+1 for num in seq1]
            else:  seq=[int(num*(N-1)) for num in seq1]
            tol = abs(mean-np.mean(seq))   
	    #print alpha, mean, np.mean(seq)       
            # check if the average of the total-degree list is close to the network mean degree (d). Tolerance=0.5 deviation from d.
            if tol>0.2:
                if mean-np.mean(seq)>0.05:alpha-=0.05
                else:alpha+=0.05
                
        return seq

##############################################################################################    

def geometric_sequence(N, mean, seqtype):
    """returns a sequence with N entries and values following a geometric distribution with mean = mean"""
    max_trial=1000
    avg = mean 
    calc_mean=0
	
    while(abs(mean-calc_mean))>0.05:
        if (mean-calc_mean)>0.05:
            avg+=0.1

        elif (mean-calc_mean)<0.05:
            avg-=0.1
        x = 1.0/avg;
        seq = [(math.log(np.random. random())/math.log(1-x)) for i in xrange(N)]
        seq=[(int(round(x))) for x in seq]
        calc_mean = sum(seq)/(1.0*len(seq)) # compute avg degree of generated degree sequence
        
    return seq
    
##############################################################################################   
