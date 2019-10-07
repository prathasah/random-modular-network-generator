#!/usr/bin/env python

"""
This module generate random modular graphs
"""
__author__ = "Pratha Sah and Shweta Bansal"
__copyright__ = "Copyright (C) 2013 Pratha Sah"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Pratha Sah"
__email__ = "ps875@georgetown.edu"


import networkx as nx 
import random as rnd
import numpy as np
#import sequence_generator as sg
import matplotlib.pyplot as plt 

#############################################################################

# Change log
# 25 July 2013: Added a function that allows variable module size according to
#               distribution modfunction defined by the user
# 10 March 2015: Fixed code. The code now generates graphs at Q=0.
#
#############################################################################

#enter the degree distribution, modularity, total network size, number of modules, and the mean degree
def generate_modular_networks(N, sfunction, modfunction, Q, m, avg_degree, verbose=False, **kwds):
    """This function generates modular random connected graph with a specified
    degree distribution and number of modules.
    Q is the desired value of modularity as defined by Newman (2004)
    n is the network size (i.e. the total number of nodes in the network)
    d is the average network degree
    sfunction is the degree distribution of the graph
    modfunction is the distribution of module size
    """
    #iterate till the module sequence is realizable
    is_valid_module_seq= False
    scale=0
    while not is_valid_module_seq:
	
        # assign nodes to modules based on module size distribution function
        mod_nodes = {}
        nc = (1.0*N)/m # average module size
        # wd_hat = average within-degree when module sizes are kept constant. This value is used to fix minimum size of modules
        wd_hat = avg_degree*(Q+(1./m))
        scale+=0.5
        mod_nodes = assign_module_membership(scale, N, nc, m, wd_hat, modfunction)
        mod_sizes = [len(mod_nodes[x]) for x in mod_nodes.keys()] 
        
        
        # calculate tolerance for module-level within-degree and average-degree
    	# 0.005 is the tolerance on Q
    	tol = 0.1 * avg_degree/(1.0* (1-sum([(num/(1.0*N))**2 for num in mod_sizes])))
        #print ("tol=="), tol
        #Qmax = 1.0 - (sum([(num/(1.0*N))**2 for num in mod_sizes]))
        #print ("Qmax=="), Qmax
        #if Q >= Qmax:
        #    raise ValueError ("Q value exceeds Qmax for the chosen number of modules. Select a lower value")    
        G = nx. Graph()
        G.add_nodes_from(range(0, N))
    
       
        # Calculate the average within-degree of the network
        wd1 = avg_degree*(Q + (sum([(num/(1.0*N))**2 for num in mod_sizes])))
        wd = round(wd1, 2)
        
        # Network generated in 5 steps
        # Step 1: Created total-degree list
        # Step 2: Create within-degree list
        # Step 3: Create between-degree list
        # Step 4: Connect outstubs (create between-edges)
        # Step 5: Connect instubs (create within-edges)
        
        connect_trial=100 # set initial connect_trial high to enter the while loop
        graph_connected = False
        outedge_graphical = False
        while connect_trial >= 10 or outedge_graphical == False or graph_connected == False:    
            connect_trial = 0
            print ("generating degree lists......")
            #assigns total-degree to each node  
            degree_list = create_total_degree_sequence (N, sfunction, avg_degree, mod_nodes, 0.01*tol, max_tries=1000, **kwds) 
            
            
            #assigns within-degree to each node                
            print ("generating indegree lists......")
            indegree_list = create_indegree_sequence(N, m, sfunction, mod_nodes, wd, degree_list, tol, **kwds)
            if len(indegree_list)==0:
                is_valid_module_seq= False
                break
            else: is_valid_module_seq= True
           
            
            #compute between-degree by formula d=wd+bd     
            outdegree_list = create_outdegree_sequence(degree_list, indegree_list) 
            # check if the outdegree (i.e between-degree) list is graphical
            outedge_graphical = is_graphical(outdegree_list, mod_nodes, m) 
            print ("outdegree list graphical?"), outedge_graphical
            if outedge_graphical == True:
                print ("connecting out nodes..............")
                #connect nodes between modules using outedge list
                connect_trial = connect_out_nodes(G, m, mod_nodes, outdegree_list,connect_trial, indegree_list) 
                
                # connect within-module edges only when between-module edges are connected.
                if len(G.edges()) > 1: 
                    print ("connecting in nodes..............")
                    #connect nodes within a module using the inedge list

                    connect_in_nodes(G, m, mod_nodes, indegree_list, outdegree_list) 
            graph_connected = nx.is_connected(G) # check if the graph is connected 

    Q1 = test_modularity_variable_mod(G, mod_nodes, verbose)
    return G
    
#############################################################################

def assign_module_membership(scale, N, nc, m, wd_hat, modfunction):
    """assigns nodal membership to module according to the specified distribution"""  
   
    mod_nodes = {}
    estimate_se = np.std(modfunction(m, nc, seqtype="modulesize"))/np.sqrt(m)
  
    min_size = int(round(wd_hat + (scale*estimate_se)))
    valid_seq = False
    while not valid_seq:
        trialseq=modfunction(m-1, nc, seqtype="modulesize") #generate m-1  number of random nos.
        seq=[max(min_size,s) for s in trialseq]
        valid_seq= (N-sum(seq)) >= min_size # so that the last module is at least size=min_size
        seq.append(N-sum(seq))
    
    init=0
    for mod in xrange(m):
        mod_nodes[mod]=[x for x in range (init, init+seq[mod])]
        init=init+seq[mod]
    
    return mod_nodes
        
#############################################################################       
    

#assign each node with degree based on user-defined distribution          
def create_total_degree_sequence (n, sfunction, avg_degree, mod_nodes, tolerance, max_tries=2000, **kwds):

        """
        Creates a total-degree sequence.Ensures that the minimum degree is 1 and
        the max degree is 1 less than the number of nodes and that the average
        degree of the sequence generated is within a tolerance of 0.05.

	`n`: number of nodes
	`sfunction`: a sequence generating function with signature (number of
			nodes, mean)
	`avg_degree`: mean degree
	`max_tries`: maximum number of tries before dying. 
	
	8Nov 2013: function now assigns total degree to the nodes module-wise, to
	ensure that mean(d(k))= mean(d). ie., mean degree of each module is equal 
	to the network mean degree
	
	"""

	
	seqlist=[]
	# this loop assumes modules are indexed sequentially from 0 to K-1
	for mod in xrange(len(mod_nodes.keys())):
		tries = 0
		max_deg = len(mod_nodes[mod]) -1
		is_valid_seq = False
        	tol = 5.0
		while tol > tolerance  or (not is_valid_seq) or (tries > max_tries):
			trialseq = sfunction(len(mod_nodes[mod]), avg_degree, seqtype="degree")
			seq = [min(max_deg, max( int(round(s)), 1 )) for s in trialseq]
			is_valid_seq = nx.is_valid_degree_sequence_havel_hakimi(seq)
			
			if not is_valid_seq and sum(seq)%2 !=0:
				x = rnd.choice(xrange(len(seq)))
				seq[x] += 1 
			is_valid_seq = nx.is_valid_degree_sequence_havel_hakimi(seq)
			# check if d_k (bar) = d(bar)
			tol = abs(avg_degree - np.mean(seq))
			tries += 1	
			if (tries > max_tries):
				raise nx.NetworkXError, \
			 	"Exceeded max (%d) attempts at a valid sequence."%max_tries
			
		seqlist.append(seq)
	deg_list = [val for sublist in seqlist for val in sublist]
	return deg_list
           
        
#############################################################################

#assign each node with within-degree based on user-defined distribution 

def create_indegree_sequence(n, m , sfunction, mod_nodes, wd,  degree_list, tolerance, **kwds):

	""" 
	Creates indegree sequence.
	Ensures that (i)the within-module degree of node is less than or equal to
	its total degree; (ii) the within-module degree sequence is graphical
	nodes; and (iii) the average within module degree of the sequence
	generated is within a tolerance of 0.05 

	'n': number of nodes
	'sfunction': a sequence generating function with signature (number of
			nodes, mean)
	 'wd' = mean within-module degree, m= total modules in the network
	 'nc'=average community size
	 'mod_nodes'= dictionary of nodal membership to communities
	'avg_degree': mean degree
	'degree_list'= total degree sequence
	"""
	is_valid_seq = False
	is_valid_indegree = False
	is_valid_module_size=False
	is_valid_outdegree = True
	tol = 5.0
	connect_trial=0
        # Return empty list if the main loop breaks
	indegree_list=[] 
	
	while ((not is_valid_seq) or (not is_valid_indegree) or (not is_valid_module_size) or (not is_valid_outdegree)):
            
            indegree_seq = sfunction(n, wd, seqtype="indegree")
            indegree_sort = list(np.sort(indegree_seq))
            degree_sort = list(np.sort(degree_list))                              
            is_valid_indegree= all([indegree_sort[i] <= degree_sort[i] for i in range(n)])==True
            mod_sizes = [len(mod_nodes[x]) for x in mod_nodes.keys()]
            is_valid_module_size = min(mod_sizes) >= max(indegree_seq)
            is_valid_seq=True
            
           
	    if not is_valid_module_size: connect_trial+=1           
            if is_valid_indegree and is_valid_seq and is_valid_module_size:
                #assign within-degree to a node such that wd(i)<=d(i)
                indegree_list = sort_inedge(indegree_sort, degree_sort, degree_list, n)
                #if network has two modules then sum(outdegree) for module 1 has to be equal to sum(outdegree) of module 2
           	if m==2: 
            		mod1, mod2 = mod_nodes.keys()
            		is_valid_outdegree = sum([(degree_list[i] -indegree_list[i]) for i in mod_nodes[mod1]]) == sum([(degree_list[i] -indegree_list[i])  for i in mod_nodes[mod2]])
            		
            	
                for module in mod_nodes:
                    seq = [indegree_list[i] for i in mod_nodes[module]]
                  
                    while (sum(seq)%2) != 0: 
                        # choose a random node in the module
                        node_add_degree = rnd.choice([i for i in mod_nodes[module]])
                        # ensure that wd<=d and wd< (module size -1)after adding a within-degree to the node
                        if indegree_list[node_add_degree] < degree_list[node_add_degree] and  indegree_list[node_add_degree] < len(mod_nodes[module]): 
                            indegree_list[node_add_degree] += 1 
                    	seq = [indegree_list[i] for i in mod_nodes[module]]
                 
                    is_valid_seq = nx.is_valid_degree_sequence_havel_hakimi(seq)
                    tol = abs(wd - (sum(seq)/(1.0*len(seq)))) #ensure that wd_k(bar) = wd(bar)
                    if (not is_valid_seq) and tol>tolerance:
                        break                     

            
            if connect_trial>10:break	    
          
        return indegree_list
   
#############################################################################
    
#assign each node with between-degree based on formula d=wd+bd   
def create_outdegree_sequence(degree_seq, indegree_seq):
    """Creates out(between-module) degree sequence"""
    outdegree_seq = [x-y for x, y in zip(degree_seq,indegree_seq)]
    return outdegree_seq
    
#############################################################################

#Connect instubs (form within-edges)
def connect_in_nodes(G, m, mod_nodes, indegree_list, outdegree_list):
    """Connects within-module stubs (or half edges) using a modified version of 
    Havel-Hakimi algorithm"""
    a = 0
   
    while a < m:
        GI = nx.Graph() # treat each module as an independent graph
        edge1 = {}
        
        #create dictionary with key=node id and value=within-degree
        for num in range (min(mod_nodes[a]), (max(mod_nodes[a])+1)):
            edge1[num] = indegree_list[num]

        # sorts within-degrees in descending order keeping the  node identity intact
        edge1_sorted = sorted(edge1.items(), key = lambda x: x[1], reverse=True) 

	# creates a tuple of (node#, within-degree) 
        nodelist, deglist = [[z[i] for z in edge1_sorted] for i in (0, 1)]
        node_edge = zip(deglist, nodelist)
        node_edge = [(x, y) for x, y in node_edge if x > 0]

        # Connection of instubs based on Havel-Hakimi algorithm.
        #Terminates when number of instubs=0
        while node_edge:
            node_edge.sort(reverse = True)
            deg1, node1 = node_edge[0] #choose node with highest within-degree
            if deg1 > 0:
                #connect the stubs of the node(i) to the next "wd(i)" nodes in the sorted list. 
                for num in range(1, 1+deg1):
                    deg2, node2 = node_edge[num]                                       
                    GI.add_edge(node1, node2)
                    node_edge = [(x-1, y) if y == node2 else (x, y) for x, y in node_edge] # reducing the degree (stubs) the next "wd(i)" nodes by 1
                    
            node_edge = node_edge[1:] #remove the node1 from the list
            node_edge = [(x, y) for x, y in node_edge if x > 0] # remove node if the within-degree (instub) hits zero
        
        randomize_graph(GI) # remove degree correlations by edge-randomization
        
        if nx.is_connected(GI) == False: # reconnect graph if disconnected
        	connect_module_graph(GI, outdegree_list)
	
	 
	#integrate the sub-graph to the main graph G.
        G.add_edges_from(GI.edges())
        a += 1
    
#############################################################################

#Connect outstubs (form between-edges)      
def connect_out_nodes(G, m, mod_nodes, outdegree_list, connect_trial, indegree_list):
    """Connects between-module stubs (or half edges) using a modified version of 
    Havel-Hakimi algorithm"""
    nbunch = G.edges()
    # additonal check: to ensure that Graph G does not have any 
    # pre-existing edges from earlier steps 
    G.remove_edges_from(nbunch) 
    is_valid_connection = False

    # maximum attemp to connect outstubs=10
    while is_valid_connection == False and connect_trial < 10:
        is_valid_connection = True     
        trial = 0
        outnodelist = []

        # creates a tuple of (node#, degree, module#)
        for num in xrange(len(G.nodes())):
            my = [key for key in mod_nodes.keys() if num in mod_nodes[key]]
            outnodelist.append((num, outdegree_list[num], my[0]))  

        #Connect outstubs using a modified version of Havel-Hakimi algorithm
        #terminate when all the outstubs are connected
        while outnodelist:
            
            if is_valid_connection == False: break
            rnd.shuffle(outnodelist)
            outnodelist = [(x, y, z) for x, y, z in outnodelist if y != 0] # removes outnodes with zero outedges
            
            # select module with the highest between-degree = hfm
            outmod_tot = [(sum([y for x, y, z in outnodelist if z == a]), a) for a in xrange(m)]
            outmod_tot.sort(reverse = True)
            hfm = outmod_tot[0][1]
	 
            # Select node (=node1) in module=hfm which has the highest between-degree
            possible_node1 = [(x, y) for x, y, z in outnodelist if z == hfm]
            possible_node1 = sorted(possible_node1, key = lambda x:x[1], reverse = True) 
            node1, deg1 = possible_node1[0] 

            
            # connect all the outstubs of node1 to random possible nodes
            # Criteria for selecting possible nodes
            # (a) Nodes cannot belong to the same module as node1
            # (b) No multi-edges allowed
            for degrees in xrange(deg1):
                if is_valid_connection == False: break
                 # criteria (a) and (b)
                node_exclude = set(G.neighbors(node1)).union(set(mod_nodes[hfm]))
                possible_node2 = [(x, y, z) for x, y, z in outnodelist if x not in node_exclude]
                # list of possible nodes that node1 can connect to.
                 
                possible_node2 = [(x, y) for x, y, z in possible_node2] 
                
                #terminate if there are no possible nodes left for node 1 to connenct to.
                if len(possible_node2) > 0:
                    is_isolates_avoided = False
                    is_valid_connection = True
                
                # prevent nodes with 0within- module edge and one between-
                #module edge to connect to one another
                #Avoids formation of disconnected graphs
                    trial=0
                    # preference to nodes with indegree>0 to connect to nodes to indegree= 0
                    if indegree_list[node1]!=0 and [indegree_list [node] for (node,deg) in possible_node2].count(0)>1: 
                    	possible_node2 = [(node, deg) for node, deg in possible_node2 if indegree_list [node]==0]
                    while not is_isolates_avoided: 
                    	    node2, deg2 = rnd.choice(list(possible_node2))
                            trial+=1
                            is_isolates_avoided = indegree_list [node1] != 0 or indegree_list[node2] != 0
                            
                            if is_isolates_avoided: trial = 0 
                            
                            # terminate if the attempt of finding possible nodes exceed 10000
                            if trial == 50000:
                           	is_valid_connection, connect_trial = remove_outedges(G, connect_trial)
                                break

                # on termination remove all the between-edges added and try again
                else:
                    is_valid_connection, connect_trial = remove_outedges(G, connect_trial)
                    break
            
                if is_valid_connection == True:
                	G.add_edge(node1, node2)
                	# reduces the degree of node1 and node 2by 1 in outnodelist
                	outnodelist = [(x, y, z) if  x != node2  else (x, y-1, z) for x, y, z in outnodelist] 
            
            # remove node if the between-degree (outstub) hits zero    	
            outnodelist = [(x, y, z) for x, y, z in outnodelist if y > 0 and x != node1] 
	  
    # remove degree correlations by edge-randomization
    if len(G.edges()) > 0:
    	randomize_graph_outedges(G, mod_nodes, indegree_list, outdegree_list)  
           
    return connect_trial                
    
#############################################################################

def remove_outedges(G, connect_trial):
    """Removes all the within-module edges of the graph"""
    edges_added = G.edges()
    G.remove_edges_from(edges_added)
    connect_trial += 1
    state = False
    
    return state, connect_trial

##################################################################3###########

def is_graphical(edgelist, mod_nodes, m):
    """Check if the between-module degree sequence is graphically realizable 
    using algorithm by Chungphaisan (1974). The algorithm allows for the 
    assumption of modules to beindividual nodes and allows multiple-edges 
    between the nodes."""
    state = True
    b = sum(edgelist)
    n = m
    s = {}
    for a in range(0, m):
        s[a+1] = sum(edgelist[min(mod_nodes[a]):(max(mod_nodes[a])+1)])
    if b%2 != 0:
        state = False
    if m==2: 
    	if s[1]!=s[2]:state=False
    else:
	    lhs=[]
	    rhs=[]
	    for j in range(1, n+1):
	    	lhs.append(sum([s[i] for i in range (1, j+1)])-(b*j*(j-1)))
	    	rhs.append(sum([min((j*b), s[i]) for i in range(j+1, n+1)]))
	   
	    if lhs > rhs:state = False
    return state

#############################################################################

#assign within-degree to a node such that wd(i)<=d(i)
def sort_inedge(inedge_sort, edge_sort, degree_list,n):
    """Assign within-module degree to nodes from the list of random-numbers 
    sequence generated with the constraint that within-module degree is less 
    than or equal to the total nodal degree."""
    indegree_list_dict = {}
    indegree_list1 = {}
    degree_list1 = {}
    nodelist = [x for x in xrange(n)]
    rnd.shuffle(nodelist) # so that node are chosen at random

    for i in range(n):
        indegree_list1[i] = inedge_sort[i]
        degree_list1[i] = edge_sort[i]
                    
    while nodelist:
        node_i = nodelist.pop() #chose a random node
        for key, value in degree_list1.items():
            #check for the rank of its total-degree in the sorted 
            #total-degree list
            if value == degree_list[node_i]: 
                #assign within-degree with the same rank (=r) in the sorted
                #within-degree list
                indegree_list_dict[node_i] = indegree_list1[key] 
                del indegree_list1[key]
                del degree_list1[key]
                break
    indegree_list = [value for key, value in indegree_list_dict.items()]
    return indegree_list
        
#############################################################################

def randomize_graph(G):
    """randomize a network using double-edged swaps.
Note: This is used to randomize only the within-edges. A separate algorithm 
(randomize_graph_outedges) is used to randomize between-edges"""

    size = G.size() # number of edges in graph
    its = 1
    for counter in range(its):
	nx.double_edge_swap(G, nswap = 5*size, max_tries = 10000*size)

#############################################################################		


def randomize_graph_outedges(G, mod_nodes, indegree_list, outdegree_list): 
    """randomize between-edges using double-edge swaps"""
    size = G.size()
    double_edge_swap_outedges(G, mod_nodes, indegree_list, outdegree_list, nswap = 5*size,max_tries = 10000*size)

#############################################################################
	
def double_edge_swap_outedges(G, mod_nodes, indegree_list, outdegree_list, nswap, max_tries):
    """Randomizes between-modul edges of the graph. This function is similiar 
    to the generic double-edge-swap technique with an additonal constraint: 
    Swaps that create within-module edges are not allowed."""
    
    if len(G) < 4: raise nx.NetworkXError("Graph has less than four nodes.")
    n = 0
    swapcount = 0
    #keys, degrees = zip(*G.degree().items()) # nodes, degree
    keys, degrees = zip(*G.degree()) # nodes, degree
    cdf = nx.utils.cumulative_distribution(degrees)  # cdf of degree
    while swapcount < nswap:
        (ui, xi) = nx.utils.discrete_sequence(2, cdistribution=cdf)
        if ui == xi :
            continue # same source, skip
        u = keys[ui] # convert index to label
        x = keys[xi]
        u_mod = [module for module in mod_nodes if u in mod_nodes[module]]
        u_mod = u_mod[0]

        if x in mod_nodes[u_mod]:
            continue # same module, skip

        # choose target uniformly from neighbors
        v = rnd.choice(list(G[u]))
        y = rnd.choice(list(G[x]))
        
        v_mod = [module for module in mod_nodes if v in mod_nodes[module]]
        v_mod = v_mod[0]
        
        if v == y or y in mod_nodes[v_mod]:
            continue # same target or same module, skip
        if (x not in G[u]) and (y not in G[v]): # don't create parallel edges
            if (indegree_list[u]+indegree_list[x]!=0)  and (indegree_list[v]+indegree_list[y]!=0):

                G.add_edge(u, x)
                G.add_edge(v, y)
                G.remove_edge(u, v)
                G.remove_edge(x, y)
                swapcount += 1
        if n >= max_tries:
            e = ('Maximum number of swap attempts (%s) exceeded '%n +
            'before desired swaps achieved (%s).'%nswap)
            raise nx.NetworkXAlgorithmError(e)
        n += 1
    return G
        
#############################################################################

def connect_module_graph(G, outdegree_list):
    """Connect disconnected modules. Note: This function cannot be used to 
    connect the entire modular graph."""
    cc_tot = list(nx.connected_components(G))  # cc returns the connected components of G as lists cc[0], cc[1], etc.  
    isolated_comp, outedge_comp, isolated_comp_count, outedge_comp_count = partition_network_components(cc_tot, outdegree_list)
    
    while isolated_comp_count > 0:   #while G is not connected, reduce number of components
        # pick a random node in the largest component cc[0] that has degree > 1
        node1 = rnd.choice(isolated_comp[0])
        # pick a node in another component whose degree >1
        node2 = rnd.choice(outedge_comp[rnd.choice([x for x in xrange(outedge_comp_count)])])
        while G.degree(node2) <= 1:
            node2 = rnd.choice(outedge_comp[rnd.choice([x for x in xrange(outedge_comp_count)])])
   
        # pick neighbors of node1 and node2
        nbr1 = rnd.choice([n for n in G.neighbors(node1)])
        nbr2 = rnd.choice([n for n in G.neighbors(node2)])

        # swap connections between node1,nbr1 with connections between node2,nbr2
        #  to attempt to connect the two components
        G.remove_edges_from([(node1, nbr1), (node2, nbr2)])
        G.add_edges_from([(node1, node2), (nbr1, nbr2)])

        cc_tot = list(nx.connected_components(G))
        isolated_comp, outedge_comp, isolated_comp_count, outedge_comp_count = partition_network_components(cc_tot, outdegree_list)		

#############################################################################

def partition_network_components(cc_tot, outdegree_list):
    """Partitions network disconnected components of a module into: 
    (a) components with no within-module edges and hence isolates, and
    (b) components with atleast one within-module edge"""
    tot_component_count = len(cc_tot)
    isolated_comp = {}
    outedge_comp = {}
    count1 = 0
    count2 = 0
    for component in xrange(tot_component_count):
            outedge_stubs = False
            for node in cc_tot[component]:
                if outdegree_list[node] > 0:
                    outedge_stubs = True
				
            if outedge_stubs is False and len(cc_tot[component]) > 1:
                isolated_comp[count1] = list(cc_tot[component])
                count1 += 1
            elif outedge_stubs is True and len(cc_tot[component]) > 1:
                outedge_comp[count2] = list(cc_tot[component])
                count2 += 1
    return isolated_comp, outedge_comp, len(isolated_comp), len(outedge_comp)
################################################################################

def adjust_indegree_seq(mod_nodes, indegree_seq):
    """checks which module does node belongs to and returns the module size"""
    mod_size=[(len(mod_nodes[s]), s) for s in mod_nodes]
    mod_size.sort()
    adjust_indegree_seq={}
    for mods in mod_size[:-1]:
        max_deg_valid=False
        total_deg_valid=False
        while not(max_deg_valid) or not(total_deg_valid):
            indegree_topop=[x for x in indegree_seq]
            rnd.shuffle(indegree_seq)
            adjust_indegree_seq[mods[1]]=[indegree_topop.pop() for nodes in mod_nodes[mods[1]]]
            max_deg_valid= max(adjust_indegree_seq[mods[1]])<mods[0]
            total_deg_valid=sum(adjust_indegree_seq[mods[1]])< mods[0]*2
        indegree_seq=[x for x in indegree_topop]   
    adjust_indegree_seq[mod_size[-1][1]]=[x for x in indegree_seq]
    indegree_seq=[]
    for mod in range(min(mod_nodes), max(mod_nodes)+1):
        for elements in adjust_indegree_seq[mod]:
            indegree_seq.append(elements)  
    return indegree_seq    
    
##############################################################################################

# 
def test_modularity_variable_mod(G, mod_nodes, verbose):
    """Computes Q value of a graph when the nodal assignment to modules is known"""
    
    mods = len(mod_nodes.keys()) # number of modules
    s= [len(mod_nodes[x]) for x in xrange (mods)]
    nodes, degrees = zip(*G.degree()) # nodes, degree
    avg_deg= np.mean(degrees)
    N=len(G.nodes())
    Q=0
    wd_bar=[]
    d_bar=[]
    

	
    for modules in range (0,mods):
        wdsum=0
        dsum=0
        wdsum=[]
        dsum=[]
        for node in mod_nodes[modules]:
            aii=[n for n in G.neighbors(node)]# total degree of a node
            a_set=set(aii) # set of all the neighbors
            mod_set=set(mod_nodes[modules]) # set of nodes in nodal module
            eii=list(a_set.intersection(mod_set)) #set of neighbors present in the same module
            wdsum.append(len(eii))
            dsum.append(len(aii))
        
        wd_bar.append(np.mean(wdsum))
        d_bar.append(np.mean(dsum))
       
    Q_list=[((wd_bar[i]*s[i])/(1.0*avg_deg*N)) - (s[i]/(1.0*N))**2 for i in xrange(mods)]
    Q= sum(Q_list)
    print ("Summary...")
    print ("averge network degree"), avg_deg
    print ("Network size"), N
    print ("List of module sizes"), s
    print ("Estimated modularity = "), Q
    if verbose:
	for modules in range (0,mods):
	    print ("Module #"+str(modules)+" has node ids ="), mod_nodes[modules]
  
    
    return Q
    
############################################################################################## 
def  generate_simple_graph(sfunction, N, avg_degree):
	"""generate a simple random graph with sfunction degree sequence"""
	
	graphical_deg_seq= False
	is_connected_graph = False
	while is_connected_graph == False:
		while graphical_deg_seq== False:
			seq = sfunction(N, avg_degree, seqtype="simple_degree")
			graphical_deg_seq= nx.is_valid_degree_sequence(seq)
		G = nx.havel_hakimi_graph(seq)
		G.remove_edges_from(G.selfloop_edges())
		   
		
		if not nx.is_connected(G):
			try: 
				connect_simple_graph(G)
				is_connected_graph = True
				randomize_graph(G)	
			except (IndexError): is_connected_graph = False
			
			
		if not nx.is_connected(G):
			try:
				connect_simple_graph(G)
				is_connected_graph = True
				
			
			except (IndexError): 
				is_connected_graph = False
				graphical_deg_seq= False
	
	return G


############################################################################################## 

def connect_simple_graph(G):
	"""check if simple graph G is disconnected and connect if necessary"""

	cc = list(nx.connected_components(G))  # cc returns the connected components of G as lists cc[0], cc[1], etc
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

		cc =  list(nx.connected_components(G))
		component_count = len(cc)  
    

#############################################################################

if __name__ == "__main__":
	
    """Main function to mimic C++ version behavior"""
    
    ### Enter the network size(N), average network degree (d), total modules in the network (m), and modularity (Q)
    N=2000
    d=10
    m=10
    Q= 0.4
    
    # specify the degree distribution of the graph. In it's current format the code can generate
    # four well known degree distribution found in biological networks - scalefree, geometric, poisson and regular distribution
    degfunction = sg.poisson_sequence

    # specify the distribution of module size. The distribution can be scalefree, geometric, poisson and regular distribution (or any aribtrary sequence)
    #in it's simplest form speicify module size tp be regular which implies that all modules are of equal size
    modfunction = sg.regular_sequence
    
    try :
    	print "Generating random modular graph with modularity(Q) = " +str(Q)
    	print "The graph has " + str(N) + " nodes, " + str(m)+ " modules, and a network mean degree of " +str(d)
    	print "Generating graph....."
   
    	G = generate_modular_networks(N, degfunction, modfunction, Q, m, d)
    	
    	print "Saving adjacency matrix under the filename adjacency_matrix.txt"
    	adjacency_matrix = nx.to_numpy_matrix(G)
	np.savetxt("adjacency_matrix.txt", adjacency_matrix, delimiter="\t")
	
    	print "writing graph in the GraphML format under the filename random_modular_graph.graphml"
   	# The graphml file can be uploaded in Gephi (http://gephi.github.io/) for graph visualization.
    	nx.write_graphml(G, "random_modular_graph.graphml")
    	
    except (IndexError, IOError):
    	print "try again"
#############################################################################

