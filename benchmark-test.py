#!/usr/bin/env python

"""
This module generates data for Figure 4, Figure S5, S6 and S7 of the manuscript.

Title: "Exploring community structure in biological networks through null modular-network graph"
Authors: Pratha Sah and Shweta Bansal

The code generates random graphs of a particular Q value and performs community detection 
on these graphs by 6 popular algorithms 

(a) Spinglass or Potts (Reichardt2006) model; 
(b) Walktrap (Pons2006) algorithm, 
(c) Infomap (Rosvall2010) algorithm, and 
(d) Label propagation model (Raghavan2007)
(e) Fast modularity optimization (Blondel2008)
(f) Fast greedy modularity optimization (Clauset2004)


Returns: 3 files each for poisson, geometric and scalefree random modular networks.
	  Each file returns the estimated Q-value by 6 community detection algorithms. 
	  For each Q-value, community detection is performed on 10 generated networks.
	  Due to stochastic nature of spingalss, infomap and label propagation models,
	  25 detection iterations are performed on each network and the average value is 
	  returned. 
	 
Requirement: 
(a)  generate-random-modular-graphs.py 
(b)  generate-random-connected-graphs-PS.py
(b)  Networkx package
(c)  Igraph package for python
"""

__author__ = "Pratha Sah and Shweta Bansal"
__copyright__ = "Copyright (C) 2013 Pratha Sah"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Pratha Sah"
__email__ = "ps875@georgetown.edu"


from igraph import *
from numpy import *
from networkx import *
from generate-random-modular-graphs.py import *
from generate-random-connected-graphs-PS.py import *

################################################################################################

# Change Log 
# 11 July 2013: Created the file
 
################################################################################################

# Enter network size (N), total modules in network (M), average network degree(avg_degree)
N=2000
M=10
avg_degree=10
graphtype1 = ["poisson", "geometric", "scalefree"]
Qrange=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

for graph1 in graphtype1:
    print ("graph="), graph1
    print ("graph="), graph1
    filename="Fig4_"+graph1+"random.txt"
    for Q in Q1:
    	print ("Q="), Q
       	for i in xrange(10):
        	c1=[]
        	spin1=[]
        	spin2=[]
        	blondel1=[]
        	infom1=[]
        	infom2=[]
        	fst1=[]
        	label1=[]
        	label2=[]
			
			# Generate random graphs
			if Q==0.0: G=generate_graph (graph1, N, avg_degree)    
        	else: G=generate_mod_networks (graph1, Q, N, M,avg_degree)
        
        	dendrogram = G.community_walktrap()
        	cl = dendrogram.as_clustering()
        	c= cl.membership
        	c1.append(G.modularity(c))
        
         	blondel=G.community_multilevel(weights=None, return_levels=False)
        	blondel1.append(G.modularity(blondel))

        	fastg=G.community_fastgreedy(weights=None)
        	fast = fastg.as_clustering()
        	fst= fast.membership
        	fst1.append(G.modularity(fst))

        	for iter in xrange(25):
            	spin= G.community_spinglass(weights=None, spins=25, parupdate=True, start_temp=2, stop_temp=0.01, cool_fact=0.99, update_rule= "config", gamma=1, lambda_=1)
            	spin1.append(G.modularity(spin))  
            
            	infom=G.community_infomap(edge_weights=None, vertex_weights=None, trials=10)
            	infom1.append(G.modularity(infom))
            
            	label=G.community_label_propagation(weights=None, initial=None, fixed=None)
            	label1.append(G.modularity(label))

        	spin2.append(mean(spin1))
        	spin2_std.append(std(spin1))
        
        	infom2.append(mean(infom1))
        	infom2_std.append(std(infom1))
        
        	label2.append(mean(label1))
        	label2_std.append(std(label1))
        	
        	f.write("Q= "+str(Q)+ '\n')
        	f.write("walktrap, "+str(c1)+'\n')
        	f.write("spinglass, "+str((spin2))+'\n')
        	f.write("spinglass std= "+ str(spin2_std)+'\n')
            f.write("blondel, "+str((blondel1))+'\n')        
        	f.write("infomap, "+str((infom2))+'\n')
        	f.write("infomap std= "+ str((infom2_std))+'\n')    
        	f.write("fast greedy, "+str((fst1))+'\n')        
        	f.write("label, "+str((label2))+'\n')
        	f.write("label std= "+ str((label2_std))+'\n')
        	
    f.close()
#################################################################################################################################
