#!/usr/bin/env python

"""
This module generates Figure 3 of the manuscript 

Title: "Exploring community structure in biological networks through null modular-network graph"
Authors: Pratha Sah and Shweta Bansal

Returns: 3 files each for poisson, geometric and scalefree random modular networks.

	 Each file returns the value of three network properties 
	 (assortativity coefficient, Clustering coefficient and Path length)
	 of random modular networks with Q range 0.1-0.8
	 
	 Each value is a list of returns from 50 iteration

Requirement: 
(a)  generate-random-modular-graphs.py 
(b)  Networkx package
"""

__author__ = "Pratha Sah and Shweta Bansal"
__copyright__ = "Copyright (C) 2013 Pratha Sah"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Pratha Sah"
__email__ = "ps875@georgetown.edu"


from networkx import *
from generate-random-modular-graphs.py  import *
from numpy import *

################################################################################################

# Change Log 
# 11 July 2013: Created the file
 
################################################################################################

# Enter the network size(n), average network degree (d), total modules in the network (m)
n=2000
d=10
m=10
Qrange=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
graphtype=["poisson", "geometric", "scalefree"]

for graph in graphtype:
        filename= "Fig3_"+graph+".txt"
        f=open(filename, "w") # open a new file for the graphtype
        for Q in Qrange:
        	assort=[]
        	cluster=[]
        	pathlength=[]
        	diam=[]
		f.write ("Q="+ str(Q)+'\n') # write Q value

        	for x in xrange(50): # 50 iteration for each Q value
			G=generate_mod_networks (graph, Q,n,m,d)
                	assort.append(degree_assortativity_coefficient(G))
                	cluster.append(average_clustering(G, nodes=None, weight=None, count_zeros=True))
                	pathlength.append(average_shortest_path_length(G, weight=None))
		f.write("assortativity= "+str(assort)+'\n')
		f.write("clustering= "+str(cluster)+'\n')
		f.write("pathlength= "+ str(pathlength)+'\n')
                

        f.close()
            
####################################################################################################		
		
		
		

