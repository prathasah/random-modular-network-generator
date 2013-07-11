#!/usr/bin/env python

"""
This module generates random connected graphs. 

"""

__author__ = "Pratha Sah and Shweta Bansal"
__copyright__ = "Copyright (C) 2013 Pratha Sah"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Pratha Sah"
__email__ = "ps875@georgetown.edu"

from networkx import *
import random as rnd
import math
from networkx.utils import pareto_sequence
from numpy import *

################################################################################################

# Change Log 
# 11 July 2013: Created the file
 
################################################################################################

# Enter network size (N), total modules in network (M), average network degree(avg_degree)
# If choosing small world, can change rewiring probability from 0.25

def generate_graph(graphtype, N, avg_degree):
	# get base graph of given graphtype
	if graphtype is "regular":
		G = random_regular_graph(avg_degree, N) # uses the networkx function
	elif graphtype is "poisson":
		G = generate_poisson_graph(N,avg_degree)  # uses the networrkx function gnm_random_graph
	elif graphtype is "geometric":
		G = generate_geometric_graph(N,avg_degree)
	elif graphtype is "scalefree":
		G= generate_scalefree_graph(N,avg_degree)
	elif graphtype is "smallworld":
		if avg_degree%2 == 0: # i.e. the average degree is even
			G = connected_watts_strogatz_graph(N, avg_degree, 0.25) # uses the networkx function
		else:
			print "Average degree cannot be odd for Small World network. Exiting. Please try again"
			return
	else:
		print "Invalid graph type chosen. Exiting. Please try again"
		return
	
	# remove self-loops
	G.remove_edges_from(G.selfloop_edges())

	# randomize graph (using double-edged swap) and make sure it's connected
        if graphtype is not "smallworld":
		if not is_connected(G):
			connect_graph(G)		
		randomize_graph(G)
		if not is_connected(G):
			connect_graph(G)
	return G

################################################################################################
def generate_poisson_graph(N, d):

	edge_list=[]
        condition=False
        while condition==False:
            condition=True
            edge_list= list(random.poisson(lam=d-1, size=N))
	    
            edge_list=[x+1 for x in edge_list] # Minimum degree allowed=1
	   
            if sum(edge_list)%2!=0:
		x=rnd.choice(range(0,N))
		edge_list[x]-=1            
            
            if abs(d-(sum(edge_list)/(1.0*len(edge_list))))>0.05:
		
                condition=False
            if is_valid_degree_sequence(edge_list, method='hh')==False:
            	
                condition=False
               
        G = havel_hakimi_graph(edge_list) 
        return G

################################################################################################
def generate_geometric_graph(N, d):

	avg_degree = d 
	calc_degree=0
	state=False
	
	while(abs(d-calc_degree))>0.1 or state==False:
		
		state=True
		if (d-calc_degree)>0.1:
                        avg_degree+=0.01
            
		elif (d-calc_degree)<0.1:
                        avg_degree-=0.01
		x = 1.0/avg_degree;
                	           
                degseq = [(math.log(random.random())/math.log(1-x)) for i in xrange(N)]
		degseq=[(int(round(x))+1) for x in degseq]
		calc_degree = sum(degseq)/(1.0*len(degseq)) # compute avg degree of generated degree sequence
                
                if is_valid_degree_sequence(degseq)==False:
                	state=False
         
	G = havel_hakimi_graph(degseq)
	return G

################################################################################################
def generate_scalefree_graph(N,d):

	alpha=N/10
	condition=False
	alpha_state=False
	while condition==False:
		condition=True	
		edgelist=[]
		edge_list1= list((random.power(alpha, size=(N))))
            	edge_list1=[abs(1-x) for x in edge_list1]
            	edge_list=[int(num*(N-2))+1 for num in edge_list1]
				
		if sum(edge_list)%2!=0:
			x=rnd.choice(range(0,N))
			edge_list[x]+=1
		
		if abs(d-(sum(edge_list)/(1.0*len(edge_list))))>0.1:
			condition=False
				
                if d-(sum(edge_list)/(1.0*len(edge_list)))>0.1:
                        alpha-=0.01
                else:
                        alpha+=0.01
        
            	if condition==True:
                        if is_valid_degree_sequence(edge_list, method='hh')==False:
                                condition=False
	G = havel_hakimi_graph(edge_list)
	return G	 

################################################################################################
def randomize_graph(G):
# randomize a network using double-edged swaps
	#print ("assortivity before swapping"),degree_assortativity_coefficient(G)
	size = G.size() # number of edges in graph
	swaps = double_edge_swap(G, nswap=5*size, max_tries= 10000*size)

################################################################################################
def connect_graph(G):
# check if G is disconnected and connect if necessary

	cc = connected_components(G)  # cc returns the connected components of G as lists cc[0], cc[1], etc.
	component_count = len(cc)
	while component_count > 1:   #while G is not connected, reduce number of components

		# pick a random node in the largest component cc[0] that has degree > 1
		node1 = rnd.choice(cc[0])
		while G.degree(node1) == 1:
			node1 = rnd.choice(cc[0])

		# pick a node in another component
		node2 = rnd.choice(cc[1])

		# pick neighbors of node1 and node2
		nbr1 = rnd.choice(G.neighbors(node1))
		nbr2 = rnd.choice(G.neighbors(node2))

		# swap connections between node1,nbr1 with connections between node2,nbr2
		#  to attempt to connect the two components
		G.remove_edges_from([(node1,nbr1),(node2,nbr2)])
		G.add_edges_from([(node1,node2),(nbr1,nbr2)])

		cc = connected_components(G)
		component_count = len(cc)

################################################################################################

if __name__ == "__main__":
	# EXAMPLE
	N=10000 # network size
	graphtype = "geometric" # network graphtype
	avg_degree = 10 # average degree
	print "Generating a simple geometric random modular graph"
    	print "Graph has 10,000 nodes, and a network mean degree of 10"
    	print "Generating graph....."
	G = generate_graph(graphtype, N, avg_degree)
	filename = "edgelist_connected_test.txt"
	write_edgelist(G, filename)
	
################################################################################################
