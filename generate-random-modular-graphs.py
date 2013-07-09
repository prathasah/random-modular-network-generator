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


from networkx import *
from random import *
from numpy import *
import matplotlib.pyplot as plt
from pylab import show
from collections import *

#################################################################################

# Change log
# 24June 2013: Corrected for the function of random number generation in scalefree distribution

################################################################################

#enter the degree distribution, modularity, total network size, number of modules, and the mean degree
def generate_mod_networks(graphtype, Q, n,m,d):
# graphtype is "regular", "poisson", "geometric", "scalefree"
# Q is the desired value of modularity as defined by Newman (2004)
# n is the network size (i.e. the total number of nodes in the network)
# d is the average network degree


    Qmax= 1.0*(m-1)/m

    if Q>=Qmax:
        raise ValueError ("Q value exceeds Qmax for the chosen number of modules. Select a lower value")    
    
    G=Graph()
    G.add_nodes_from(range(0,n))
    graph=[]
    n_nodes=[]
    # mod nodes stores the the nodes presents in a community in a dictionary
    mod_nodes={}
    init_nodes={}
    e1=[]
    e={}

    #nodes of each module=nc, Number of nodes in modules (i.e nc) is constant
    nc=(n/m)

    # Calculate the average within-degree of the network
    wd1= (((1.0*d/m)*((Q*m+1))))
    wd=round(wd1,2)
 
    #Assign nodes to modules.
    count=0
    for a in range (0,m):
        mod_nodes[a]=range(count, (count+nc))
        count=count+nc
        
    # Network generated in 5 steps
    # Step 1: Created total-degree list
    # Step 2: Create within-degree list
    # Step 3: Create between-degree list
    # Step 4: Connect outstubs (create between-edges)
    # Step 5: Connect instubs (create within-edges)
    connect_trial=0
    graph_connected=False
    
    while connect_trial==0 or connect_trial==10 or outedge_graphical==False or graph_connected==False:
        
        connect_trial=1
        edges_intact=True 
        edge_list=assign_network_degree(graphtype, mod_nodes, m,nc,d) #assigns total-degree to each node
        print ("sum of  edge list"), sum(edge_list)
            
        inedge_list=assign_indegrees(graphtype, mod_nodes,m,nc,wd,d,edge_list)  #assigns within-degree to each node
        print ("sum inedge list"), sum(inedge_list)
               
        outedge_list=assign_outdegrees(m,nc,edge_list, inedge_list) #compute between-degree by formula d=wd+bd
        print ("sum outedge list"), sum(outedge_list)
        
        outedge_graphical=is_graphical(outedge_list,mod_nodes, m,nc) # check if the outdegree (i.e between-degree) list is graphical
        
        if outedge_graphical==True:
            
            print ("connecting out nodes")
            connect_trial=connect_out_nodes(G,m,nc, edge_list, mod_nodes, outedge_list,connect_trial,inedge_list,graphtype) #connect nodes between modules using outedge list
            
	    
            if len(G.edges())>1:
                print ("connecting in nodes")
                connect_in_nodes(G,m,edge_list, mod_nodes, inedge_list,outedge_list) #connect nodes within a module using the inedge list

        graph_connected=is_connected(G) # check if the graph is connected
        
    return G

#########################################################################################

#assign each node with degree based on user-defined distribution          
def assign_network_degree(graphtype, mod_nodes, m,nc,d): 
    
    if graphtype is "regular":
        edge_list=[d]*(m*nc)
        return edge_list

    if graphtype is "poisson":
        edge_list=[]
        condition=False
        while condition==False:
            condition=True
            edge_list= list(random.poisson(lam=d, size=(m*nc)))
            edge_list=[1 if x==0 else x for x in edge_list] # replaces all occurnces of zero degree with degree 1

            #check if the sum of total-degree in each module is even. Add an edge to a random node in the module if sum is odd
            for a in range(0,m):
                if sum(edge_list[min(mod_nodes[a]):(max(mod_nodes[a])+1)])%2!=0: 
                    x=choice(mod_nodes[a])
                    edge_list[x]=edge_list[x]+1 
			            
            # check if the average of the total-degree list is close to the network mean degree (d). Tolerance=0.5 deviation from d.
            if abs(d-(sum(edge_list)/(1.0*len(edge_list))))>0.05:
                condition=False
      
        return edge_list

    if graphtype is "geometric":

        edge_list=geometric_seq(d,m,nc, seq="tot_degree")
        #check if the sum of total-degree in each module is even. Add an edge to a random node in the module if sum is odd
        for a in range(0,m):
            if sum(edge_list[min(mod_nodes[a]):(max(mod_nodes[a])+1)])%2!=0: 
                x=choice(mod_nodes[a])
                edge_list[x]=edge_list[x]+1

        return edge_list

    if graphtype is "scalefree":
       
        alpha=(m*nc)/10
        condition=False
        while condition==False:
            condition=True
            edge_list=[]
            edge_list1= list((random.power(alpha, size=(m*nc))))
            edge_list1=[abs(1-x) for x in edge_list1]
            edge_list=[int(num*(m*nc-2))+1 for num in edge_list1]
            
            #check if the sum of total-degree in each module is even. Add an edge to a random node in the module if sum is odd
            for a in range(0,m):
                if sum(edge_list[min(mod_nodes[a]):(max(mod_nodes[a])+1)])%2!=0: 
                    x=choice(mod_nodes[a])
                    edge_list[x]=edge_list[x]+1
           
            # check if the average of the total-degree list is close to the network mean degree (d). Tolerance=0.5 deviation from d.
            if abs(d-(sum(edge_list)/(1.0*len(edge_list))))>0.05:
                condition=False
                if d-(sum(edge_list)/(1.0*len(edge_list)))>0.05:
                    alpha-=0.01
                else:
                    alpha+=0.01
        
        return edge_list
        
#############################################################################################

#assign each node with within-degree based on user-defined distribution 
def assign_indegrees(graphtype, mod_nodes, m,nc,wd,d,edge_list):
    condition=False
    while condition==False:
        inedge_list=[]
        condition=True

        if graphtype=="regular":

            if (2.0*wd)-d>0:
                inedge_list=list(random.random_integers((wd-(d-wd)),high=(wd+(d-wd)),size=m*nc))
            else:
                inedge_list=list(random.random_integers(0,high=(2*wd),size=m*nc))
        
        if graphtype=="poisson":
            state=False
            while state==False:
                state=True
                inedge_list2=random.poisson(lam=wd,size=m*nc)
                inedge_sort=list(sort(inedge_list2))
                edge_sort=list(sort(edge_list))
                diff=[x-y for x,y in zip(edge_sort,inedge_sort)]
                
                for i in range(m*nc):
                    if inedge_sort[i]>edge_sort[i]:
                        state=False
            #assign within-degree to a node such that wd(i)<=d(i)
            inedge_list=sort_inedge(inedge_sort,edge_sort,edge_list, m,nc)

        if graphtype=="geometric":
            state=False
            while state==False:
                state=True
                inedge_list2=geometric_seq(wd,m,nc, seq="in_degree")
                
                inedge_sort=list(sort(inedge_list2))
                
                edge_sort=list(sort(edge_list))
                for i in range(m*nc):
                    if inedge_sort[i]>edge_sort[i]:
                        state=False
                                                
            inedge_list=sort_inedge(inedge_sort,edge_sort,edge_list, m,nc)

        if graphtype=="scalefree":
            state=False
            alpha=nc/5
            while state==False:
                state=True
                inedge_list=[]
                inedge_list1= list((random.power(alpha, size=(m*nc)))) #gives a random probability value from 0 to 1 that follows power law
                inedge_list1=[abs(1-x) for x in inedge_list1]
                inedge_list2=[int(num*(nc-1)) for num in inedge_list1] #converting 0-1 scale to 0 to (nc-1) scale
                
                if abs(wd-(sum(inedge_list2)/(1.0*len(inedge_list2))))>0.1:
                    state=False
                    if wd-(sum(inedge_list2)/(1.0*len(inedge_list2)))>0.1:
                        alpha-=0.01
                    else:
            	        alpha+=0.01
                
                if state!=False:
                    inedge_sort=list(sort(inedge_list2))
                    edge_sort=list(sort(edge_list))
                
                    for i in range(m*nc):
                        if inedge_sort[i]>edge_sort[i]:
                            state=False
                            
            inedge_list=sort_inedge(inedge_sort,edge_sort,edge_list, m,nc)
            
        #check if the sum of within-degree in a module is even. Add an indegree to a random node if sum is odd.
        for a in xrange(m):
            if sum(inedge_list[min(mod_nodes[a]):(max(mod_nodes[a])+1)])%2!=0: 
                stat=False
                while stat==False:
                    node_add_edge=randint((min(mod_nodes[a])),(max(mod_nodes[a])))
                    if inedge_list[node_add_edge]<edge_list[node_add_edge]: # ensure that wd<=d after adding a within-degree to the node
                        stat=True
                        inedge_list[node_add_edge]+=1 

            #checks if module within-degree sequence is graphical 
            if is_valid_degree_sequence (inedge_list[min(mod_nodes[a]):(max(mod_nodes[a])+1)])==False: 
                condition=False

        # check if the average of the within-degree list is close to the network mean within-degree (wd). Tolerance=0.5 deviation from wd.
        if abs(wd-(sum(inedge_list)/(1.0*len(inedge_list))))>0.05:
            condition=False
            
    return inedge_list
   
#####################################################################################################
    
#assign each node with between-degree based on formula d=wd+bd   
def assign_outdegrees(m,nc,edge_list, inedge_list):
    outedge_list=[x-y for x,y in zip(edge_list,inedge_list)]
    return outedge_list
    
#####################################################################################################

#Connect instubs (form within-edges)
def connect_in_nodes(G,m,edge_list, mod_nodes, inedge_list,outedge_list):
    a=0
    count_inedge=0  # for code check. count the number of outedge connections made
    while a<m:
        GI=Graph() # treat each module as an independent graph
        trial=0
        edge1={}
        
        #create dictionary with key=node id and value=within-degree
        for num in range (min(mod_nodes[a]),(max(mod_nodes[a])+1)):
            edge1[num]=inedge_list[num]

        # sorts within-degrees in descending order keeping the  node identity intact
        edge1_sorted=sorted(edge1.items(), key=lambda x: x[1], reverse=True) 

	# creates a tuple of (node#, within-degree) 
        nodelist,deglist= [[z[i] for z in edge1_sorted] for i in (0,1)]
        node_edge=zip(deglist,nodelist)
        node_edge=[(x,y) for x,y in node_edge if x>0]

        # Connection of instubs based on Havel-Hakimi algorithm.
        #Terminates when number of instubs=0
        while node_edge:
            node_edge.sort(reverse=True)
            deg1,node1=node_edge[0] #choose node with highest within-degree
            if deg1>0:
                #connect the stubs of the node(i) to the next "wd(i)" nodes in the sorted list. 
                for num in range(1, 1+deg1):
                    deg2,node2=node_edge[num]                                       
                    GI.add_edge(node1,node2)
                    count_inedge+=1 
                    node_edge=[(x-1,y) if y==node2 else (x,y) for x,y in node_edge] # reducing the degree (stubs) the next "wd(i)" nodes by 1
                    
            node_edge=node_edge[1:] #remove the node1 from the list
            node_edge=[(x,y) for x,y in node_edge if x>0] # remove node if the within-degree (instub) hits zero
        
        randomize_graph(GI) # remove degree correlations by edge-randomization
        
        if is_connected(GI)==False: # reconnect graph if disconnected
        	connect_module_graph(GI,outedge_list)
	
	#integrate the sub-graph to the main graph G.
        G.add_edges_from(GI.edges())
        a+=1
    print ("number of inedge connections="), count_inedge
    
#####################################################################################################

#Connect instubs (form within-edges)      
def connect_out_nodes(G,m,nc, edge_list, mod_nodes, outedge_list,connect_trial,inedge_list,graphtype):
    nbunch=G.edges()
    G.remove_edges_from(nbunch) # additonal check: to ensure that Graph G does not have any edges from previous steps
    state=False

    # maximum attemp to connect outstubs=10
    while state==False and connect_trial<10:
        state=True     
        trial=0
        edge1=[]
        outnodelist=[]
        edge1={}

        # creates a tuple of (node#, degree, module#)
        for num in xrange(m*nc):
            my=[key for key in mod_nodes.keys() if num in mod_nodes[key]]
            outnodelist.append((num,outedge_list[num], my[0]))  

        #Connect outstubs using a modified version of Havel-Hakimi algorithm
        #terminate when all the outstubs are connected
        while outnodelist:
            if state==False:
                break
            shuffle(outnodelist)
            outnodelist=[(x,y,z) for x,y,z in outnodelist if y!=0] # removes outnodes with zero outedges
            hfm=0
            hfm_tot=0
            
            # select module with the highest between-degree = hfm
            for a in xrange(m):
            	outmod_tot= sum([y for x,y,z in outnodelist if z==a])
            	if outmod_tot>hfm_tot:
            		hfm_tot=outmod_tot
            		hfm=a

            # Select node (=node1) in module=hfm which has the highest between-degree
            possible_node1=[]
            possible_node1=[(x,y) for x,y,z in outnodelist if z==hfm]
            possible_node1=sorted(possible_node1, key=lambda x:x[1], reverse=True) 
            node1,deg1=possible_node1[0] 
            
            # connect all the outstubs of node1 to random possible nodes
            # Criteria for selecting possible nodes
            # (a) Nodes cannot belong to the same module as node1
            # (b) No multi-edges allowed
            for degrees in xrange(deg1):
                if state==False: break
                node_exclude=set(G.neighbors(node1)).union(set(mod_nodes[hfm])) # criteria (a) and (b)
                possible_node2=[(x,y,z) for x,y,z in outnodelist if x not in node_exclude]
                possible_node2=[(x,y) for x,y,z in possible_node2] # is the list of possible nodes that node1 can connect to.
                #terminate if there are no possible nodes left for node 1 to connenct to.
                if len(possible_node2)>0:
                    cd=False
                    state=True
                # prevent nodes with 0 inedge and 1 outedge to connect to one another -->Avoids formation of disconnected graphs
                    while cd==False: 
                            cd=True
                            trial==0
                            node2,deg2= choice(list(possible_node2))
                            if inedge_list [node1]==0 and inedge_list[node2]==0:
                                cd=False
                                trial+=1
                            # terminate if the attempt of finding possible nodes exceed 10000
                            if trial==10000:
                                edges_added=G.edges()
                                G.remove_edges_from(edges_added)
                                connect_trial+=1
                                state=False
                                break

                # on termination remove all the between-edges added and try again
                else:
                    edges_added=G.edges()
                    G.remove_edges_from(edges_added)
                    connect_trial+=1
                    state=False
                    break
            
                if state==True:
                	G.add_edge(node1,node2)
                	outnodelist=[(x,y,z) if  x!=node2  else (x,y-1,z) for x,y,z in outnodelist] # reduces the degree of node1 and node 2by 1 in outnodelist
                	
            outnodelist=[(x,y,z) for x,y,z in outnodelist if y>0 and x!=node1] # remove node if the between-degree (outstub) hits zero

    print ("checking from the graph outconnections are"), len(G.edges())

    # remove degree correlations by edge-randomization
    if len(G.edges())>0:
    	randomize_graph_outedges(G, mod_nodes, inedge_list, outedge_list)  
           
    return connect_trial                
    
##############################################################################################

# check if the degree sequence is graphically realizable using algorithm by Chungphaisan (1974)
def is_graphical(edgelist,mod_nodes, m,nc):
    state=True
    b=sum(edgelist)
    n=m
    s={}
    for a in range(0,m):
        s[a+1]=sum(edgelist[min(mod_nodes[a]):(max(mod_nodes[a])+1)])
    if b%2!=0:
        state=False
    for j in range(1,n):
        lhs=sum((s[i]-(b*j*(j-1))) for i in range (1,j+1))
        rhs=sum((min((j*b),s[i])) for i in range(j+1, n+1))
        if lhs>rhs:
            state=False
    return state

##############################################################################################    
   
# generate geometric degree sequence
def geometric_seq(d,m,nc, seq):
    N=m*nc
    max_trial=1000
    avg_degree = d +0.21
    calc_degree=0
	
    while(abs(d-calc_degree))>0.1:
        if (d-calc_degree)>0.1:
            avg_degree+=0.1

        elif (d-calc_degree)<0.1:
            avg_degree-=0.1
        x = 1.0/avg_degree;
        degseq = [(math.log(random.random())/math.log(1-x)) for i in xrange(N)]
        degseq=[(int(round(x))) for x in degseq]
        if seq=="tot_degree":
                degseq=[x+1 for x in degseq]
        calc_degree = sum(degseq)/(1.0*len(degseq)) # compute avg degree of generated degree sequence
        
    return degseq

##############################################################################################  

#assign within-degree to a node such that wd(i)<=d(i)
def sort_inedge(inedge_sort,edge_sort,edge_list, m,nc):
    inedge_list=[]
    inedge_list_dict={}
    inedge_list1={}
    edge_list1={}
    nodelist=[x for x in range(m*nc)]
    shuffle(nodelist) # so that node are chosen at random

    for i in range(m*nc):
        inedge_list1[i]=inedge_sort[i]
        edge_list1[i]=edge_sort[i]
                    
    iter=0
    while nodelist:
        iter=nodelist.pop() #chose a random node
        for key,value in edge_list1.items():
            if value==edge_list[iter]: #check for the rank (=r)of its total-degree in the sorted total-degree list
                inedge_list_dict[iter]=inedge_list1[key] # assign within-degree with the same rank (=r) in the sorted within-degree list
                iter+=1
                del inedge_list1[key]
                del edge_list1[key]
                break
    inedge_list=[value for key,value in inedge_list_dict.items()]
    return inedge_list
        
##############################################################################################

# Computes Q value of a graph when the nodal assignment to modules is known
def test_modularity (G,n,m):

    dict = {}
    nc=n/m
    mod_nodes={}
    count=0
    for a in range (0,m):
        mod_nodes[a]=range(count, (count+nc))
        count=count+nc
    Q=0
    for modules in range (0,m):
        asum=0
        esum=0
        eii=[]
        for node in mod_nodes[modules]:
            aii=(G.neighbors(node))# total degree of a node
            a_set=set(aii)
            mod_set=set(mod_nodes[modules])
            eii=list(a_set.intersection(mod_set))
            asum=asum+ len(aii)
            esum=esum+len(eii)
       
        denom=2.0*len(G.edges())
        Q_calc=((esum/denom)-((asum/denom)**2))
        Q=Q+Q_calc

    return Q
    
##############################################################################################

# randomize a network using double-edged swaps.
#Note: This is used to randomize only the within-edges. A separate algorithm (randomize_graph_outedges) is used to randomize between-edges
def randomize_graph(G):

	size = G.size() # number of edges in graph
	its = 1
	for counter in range(its):
		swaps = double_edge_swap(G,nswap=5*size, max_tries= 10000*size)

##############################################################################################		

# randomize between-edges using double-edge swaps.
def randomize_graph_outedges(G, mod_nodes, inedge_list, outedge_list): 
	
	size = G.size()
	swaps = double_edge_swap_outedges(G,mod_nodes, inedge_list, outedge_list, nswap=5*size, max_tries= 10000*size)
	
# similiar to the generic double-edge-swap technique with additonal constraint.
# Constraint: Swaps that create within-edges instead of between-edges are not allowed.
def double_edge_swap_outedges(G,mod_nodes, inedge_list, outedge_list, nswap, max_tries):
    if len(G) < 4:
        raise nx.NetworkXError("Graph has less than four nodes.")
    n=0
    swapcount=0
    keys,degrees=zip(*G.degree().items()) # keys, degree
    cdf=utils.cumulative_distribution(degrees)  # cdf of degree
    while swapcount < nswap:
        (ui,xi)=utils.discrete_sequence(2,cdistribution=cdf)
        if ui==xi :
            continue # same source, skip
        u=keys[ui] # convert index to label
        x=keys[xi]
        u_mod=[module for module in mod_nodes if u in mod_nodes[module]]
        u_mod=u_mod[0]

        if x in mod_nodes[u_mod]:
            continue # same module, skip

        # choose target uniformly from neighbors
        v=choice(list(G[u]))
        y=choice(list(G[x]))
        
        v_mod=[module for module in mod_nodes if v in mod_nodes[module]]
        v_mod=v_mod[0]
        
        if v==y or y in mod_nodes[v_mod]:
            continue # same target or same module, skip
        if (x not in G[u]) and (y not in G[v]): # don't create parallel edges
            if (inedge_list[u]+inedge_list[x]!=0)  and (inedge_list[v]+inedge_list[y]!=0):

                G.add_edge(u,x)
                G.add_edge(v,y)
                G.remove_edge(u,v)
                G.remove_edge(x,y)
                swapcount+=1
        if n >= max_tries:
            e=('Maximum number of swap attempts (%s) exceeded '%n +
            'before desired swaps achieved (%s).'%nswap)
            raise nx.NetworkXAlgorithmError(e)
        n+=1
    return G
        
##############################################################################################

# Connect disconnected modules.
# Note: The module cannot be used to connect the entire graph.
def connect_module_graph(G,outedge_list):
    cc_tot = connected_components(G)  # cc returns the connected components of G as lists cc[0], cc[1], etc.
    tot_component_count = len(cc_tot)
    cc={}
    count1=0
    count2=0
    cc2={}
	
    for x in xrange(tot_component_count):
        outedge_stubs=False
        for node in cc_tot[x]:
            if outedge_list[node]>0:
                outedge_stubs=True
            if outedge_stubs==False and len(cc_tot[x])>1:
                cc[count1]=cc_tot[x]
                count1+=1
            elif outedge_stubs==True and len(cc_tot[x])>1:
                cc2[count2]=cc_tot[x]
                count2+=1
    component_count = len(cc)
    component_cc2=len(cc2)

    while component_count > 0:   #while G is not connected, reduce number of components
        # pick a random node in the largest component cc[0] that has degree > 1
        node1 = choice(cc[0])
        # pick a node in another component whose degree >1
        node2 = choice(cc2[choice([x for x in xrange(component_cc2)])])
        while G.degree(node2)<=1:
            node2 = choice(cc2[choice([x for x in xrange(component_cc2)])])
   
        # pick neighbors of node1 and node2
        nbr1 = choice(G.neighbors(node1))
        nbr2 = choice(G.neighbors(node2))

        # swap connections between node1,nbr1 with connections between node2,nbr2
        #  to attempt to connect the two components
        G.remove_edges_from([(node1,nbr1),(node2,nbr2)])
        G.add_edges_from([(node1,node2),(nbr1,nbr2)])

        cc_tot = connected_components(G)
        tot_component_count = len(cc_tot)
        cc={}
        cc2={}
        count1=0
        count2=0
        
        for x in xrange(tot_component_count):
            outedge_stubs=False
            for node in cc_tot[x]:
                if outedge_list[node]>0:
                    outedge_stubs=True
				
            if outedge_stubs==False and len(cc_tot[x])>1:
                cc[count1]=cc_tot[x]
                count1+=1
            elif outedge_stubs==True and len(cc_tot[x])>1:
                cc2[count2]=cc_tot[x]
                count2+=1
		
        component_count=len(cc)
        component_cc2=len(cc2)		
		
##############################################################################################
