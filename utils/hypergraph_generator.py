''' This  file contains algorithms which ganerate random hypergraphs for given parameters. 

    H=(V,E) pair is a hypergraph

    The hypergraph implementation is described in hypergraphs.py
    
    V: vertices
    E: hyperedges (subsets of V)
    
    Hypergraph types:
    
    GIRH: Geometric Inhomogenious Random Hypergraph
    
    V: is a finite ranom set of an d-dimensional points in the d-dimensioanl cube
    
    E: m hyperedge generation operates as follows:
        Every node has a w_v weight. 
        We generate hyperedges step by step. 
        In the step i, we pick a random point in the n-dimensioanl cube.
        This point will be the hyperedge i central point c_i. After that,
        we generate the hyperedge weight w_e according to some hyperedge size distribution.
        And we chose the nodes which are contianed by the hyperede with probability proportinal:
        
        w_e*w_v/(W*||c_i-v||_2)^(d*alpha),
        
        where alpha is a parameter given as input. W is the sum of the weights.
'''

import numpy as np
import random as rnd
import math
import powerlaw
from scipy.stats import poisson
from scipy.stats import pareto
from scipy.stats import gamma
from scipy.spatial import distance
import bisect
import hypernetx as hnx
from utils.GIRH import GIRH
from collections import defaultdict


def get_hypergraph(type, args):
    if type=='girh':
        H=get_girh(args['n'],args['m'],args['alpha'], args['tau'], args)
        return H
    if type=='flower':
        H=get_flower(args)
        return H
    if type=='erdos-renyi':
        H=get_erdosrenyi(args)
        return H
    if type=='regular':
        H=get_regular(args)
        return H
    if type=='barabasi':
        H=get_preferential_attachment(args)
        return H
    if type=='toy_country':
        H, edge_types=get_toy_country(args)
        return H,edge_types
    if type=='configuration':
        H=get_config(args)
        return H

def get_toy_country(args):
    list_of_edges = {}
    edge_types=defaultdict(lambda:'f')
    n=args['n']
    edge_ind = 1
    nodelist=np.arange(start=1,stop=n+1)
    
    #Generate household edges
    if args['distribution'] == 'uniform':
        for i in range(int(n/args['house size'])):
            if (i+1)*args['house size'] > n:
                list_of_edges[str(edge_ind)] = np.copy(nodelist[(i*args['house size']):n])
            else:
                list_of_edges[str(edge_ind)] = np.copy(nodelist[(i*args['house size']):((i+1)*args['house size'])])
            edge_ind += 1
    
    #Generate centrum edge
    if 'centrum size' not in args.keys():
        centrum_size = int(0.05 * n)
    else:
        centrum_size = args['centrum size']
    centrum_size=int(centrum_size/args['house size'])*args['house size'] 
    centrum = nodelist[0:centrum_size]
    list_of_edges[str(edge_ind)] = np.copy(centrum)
    edge_types[str(edge_ind)]='centrum'
    centrum_ind=edge_ind
    edge_ind += 1
    cities=[0]
    
    
    #Generate cities
    if 'city size' not in args.keys():
        city_size = 0.01 * n
    else:
        city_size = args['city size']
    int(city_size/args['house size'])*args['house size']
    for i in range(args['city num']):
        list_of_edges[str(edge_ind)] = np.copy(nodelist[centrum_size+i*city_size:centrum_size+(i+1)*city_size])
        edge_types[str(edge_ind)]='city'
        cities.append(i+1)
        edge_ind += 1
        
        
    #Generate settlements
    begin = centrum_size+args['city num']*city_size
    settlement_size=20
    for i in range(int((n-begin)/settlement_size)):
        if begin+(i+1)*settlement_size > n:
            list_of_edges[str(edge_ind)] = np.copy(nodelist[begin+(i*settlement_size):n])
        else:
            list_of_edges[str(edge_ind)] = np.copy(nodelist[begin+(i*settlement_size):begin+((i+1)*settlement_size)])
        edge_types[str(edge_ind)]='s'
        edge_ind += 1
        cities.append(i+1+args['city num'])
        
        
    #Generate workplaces
    #Decide where the persons work
    list_of_edges[str(edge_ind)] = np.copy(nodelist)
    edge_types[str(edge_ind)]='country'
    edge_ind += 1
    
    prob_work = args['workplace prob'][0]
    workplaces=np.random.choice(cities,size = centrum_size, replace=True, p = prob_work)
    for i in range(args['city num']):
        prob_work = args['workplace prob'][i+1]
        workplace=np.random.choice(cities,size= city_size, replace=True, p = prob_work)
        workplaces=np.concatenate((workplaces, workplace), axis=None)
    
    for i in range(int((n-begin)/settlement_size)):
        prob_work = args['workplace prob'][i+1+args['city num']]
        workplace=np.random.choice(cities,size= settlement_size, replace=True, p = prob_work)
        workplaces=np.concatenate((workplaces, workplace), axis=None)
    
    for city in cities:
        workers=np.copy(nodelist[workplaces==city])
        np.random.shuffle(workers)
        
        if args['distribution'] == 'uniform':
            for i in range(int(len(workers)/args['workplace size'])):
                if (i+1)*args['workplace size'] > len(workers):
                    list_of_edges[str(edge_ind)] = np.copy(workers[(i*args['workplace size']):len(workers)])
                else:
                    list_of_edges[str(edge_ind)] = np.copy(workers[(i*args['workplace size']):((i+1)*args['workplace size'])])
                edge_types[str(edge_ind)]='w'
                edge_ind += 1
                
    
    
    
    H=hnx.Hypergraph(list_of_edges)
    return H,edge_types
    
    
def get_erdosrenyi(args):
    list_of_edges = {}
    n=args['n']
    i = 0
    nodes=np.arange(start=0,stop=n)
    while i < args['m']:
        nodelist=[]
        j = 0
        while j < args['size']:
            node=rnd.choice(nodes)
            if node not in nodelist:    
                nodelist.append(node)
                j += 1
        if nodelist not in list_of_edges.values():
            list_of_edges[str(i)]=nodelist
            i += 1
    return hnx.Hypergraph(list_of_edges)
    
    
def get_girh(n,m, alpha, tau, args):
    G=GIRH(n,m,alpha,tau,
      expected_node_weight=args['expected_node_weight'],
      expected_edge_weight=args['expected_edge_weight'], 
      C_1=1,
      C_2=1,
      on_torus=False,
      dimension=2)
    return G.hypergraph
    
def get_flower(args):
    k=args["edge size"]
    s=args["number of levels"]
    scenes={str(0) : tuple(np.arange(k))}
    m=int(k*(((k-1)**(s-1))-1)/(k-2))
    n=int(k*(((k-1)**(s))-1)/(k-2))
    v=list(np.arange(n))
    for i in range(m):
        scenes[str(i+1)]=tuple([i]+v[(i+1)*(k-1)+1:(i+2)*(k-1)+1])
    H=hnx.Hypergraph(scenes)
    
    return H  
    
def get_regular(args):
        """
        Generate random family hyperedges from hyperedge distribution.
        parameter:
        -------
        args: dict
        parameter of the hypergraph
        
        """
        n=args['n']
        list_of_edges = {}
        edge_ind = 0
        nodelist=np.arange(start=1,stop=n+1)
        if args['distribution']=='uniform':
            for k in range(args['d']):
                np.random.shuffle(nodelist)
                for i in range(int(n/args['size'])):
                    if (i+1)*args['size'] > n:
                        list_of_edges[str(edge_ind)]=np.copy(nodelist[(i*args['size']):n])
                    else:
                        list_of_edges[str(edge_ind)]=np.copy(nodelist[(i*args['size']):((i+1)*args['size'])])
                    edge_ind+=1
        return hnx.Hypergraph(list_of_edges)

def get_preferential_attachment(args):
    n=args['n']
    m = args['m']
    list_of_edges = {}
    edge_ind = 0
    nodelist=np.arange(n)
    k=0
    while k < 2:
        k =sample(args)
    list_of_edges[str(edge_ind)]=nodelist[0:k]
    edge_ind +=1
    degree_vec = np.ones(k)
    node_indexes = np.arange(k,dtype=int)
    prob_vec = degree_vec/sum(degree_vec)
    t = k
    for node in nodelist[k:]:
        for i in range(m):
            
            #sample the size of the new hyperedge
            
            size=0
            while size<1:
                size = sample(args)
            if t < size:
                size=t
            
            #sample the nodes for the new hyperedge proportional to their degrees
            
            new_hyperedge = np.random.choice(node_indexes,size= size-1, replace=False, p= prob_vec)
            
            #update the degree vector
            degree_vec[new_hyperedge] +=1
            prob_vec = degree_vec/sum(degree_vec)
            #add the hyperedge
            list_of_edges[str(edge_ind)] = np.append(new_hyperedge, node)
            edge_ind +=1
        node_indexes=np.append(node_indexes,t)
        degree_vec=np.append(degree_vec,m)
        prob_vec = degree_vec/sum(degree_vec)
        t += 1
    return hnx.Hypergraph(list_of_edges)
    
def get_config(args):
    stubs=[]
    n=args['n']
    node_ind=0
    edge_ind=0
    for node_degree in args['degree sequence']:
        for i in range(node_degree):
            stubs.append(node_ind)
        node_ind += 1
    rnd.shuffle(stubs)
    list_of_edges = generate_nodelists(edge_ind, stubs, args)
    return hnx.Hypergraph(list_of_edges)
    
    
def sample(distribution_data):
        '''
        Returns one sample element from the distribution defined in the distribution_data.

        The distribution_data can be given as binomial, poisson, uniform, 
        or a finite distribution vector.
        '''
        if distribution_data['distribution']=='binom':
            p = 1 - (distribution_data['size variance'] / distribution_data['size mean'])
            n = distribution_data['size mean'] / p
            return np.random.binomial(n, p, 1)[0]
        
        elif distribution_data['distribution']=='poisson':
            mu = distribution_data['size mean']
            return poisson.rvs(mu)

        elif distribution_data['distribution']=='gamma':
            a = distribution_data['a']
            scale = distribution_data['scale']
            loc = distribution_data['loc']
            return int(gamma.rvs(a,scale=scale,loc=loc))
        
        elif distribution_data['distribution']=='pareto':
            b = distribution_data['b']
            scale = distribution_data['scale']
            loc = distribution_data['loc']
            return int(pareto.rvs(b,scale=scale,loc=loc))

        elif distribution_data['distribution']=='uniform':
            return distribution_data['size']

        else:
            r = rnd.random()
            return bisect.bisect(distribution_data['distribution'],r)-1
        
def generate_nodelists(edge_ind, stack, hyperedge_data):
        """
        Generate random family hyperedges from hyperedge distribution.
        parameter:
        -------
        _list: dict
        dictionary of the nodes
        
        """
        list_of_edges = {}
        if hyperedge_data['distribution']=='uniform':
            for i in range(int(len(stack)/hyperedge_data['size'])):
                if (i+1)*hyperedge_data['size'] > len(stack):
                    list_of_edges[str(edge_ind)]=stack[(i*hyperedge_data['size']):len(stack)]
                    edge_ind+=1
                else:
                    list_of_edges[str(edge_ind)]=stack[(i*hyperedge_data['size']):(((i+1)*hyperedge_data['size']))]
                    edge_ind+=1
            
        if hyperedge_data['distribution']=='poisson':
            mu = hyperedge_data['size mean']
            
            while stack != []:
                size = poisson.rvs(mu)
                if size>0:
                    family = []
                    if len(stack) < size:
                        family = list(stack)
                        for i in family:
                            stack.remove(i)
                    else:
                        i=0
                        while i < size:
                            member = rnd.choice(stack)
                            if member not in family:
                                stack.remove(member)
                                family.append(member)
                                i+=1
                    list_of_edges[str(edge_ind)]=family
                    edge_ind+=1
                
        elif hyperedge_data['distribution']=='binom':
            p = 1 - (hyperedge_data['size variance'] / hyperedge_data['size mean'])
            n = hyperedge_data['size mean'] / p
            
            while stack != []:
                size = np.random.binomial(n, p, 1)[0]
                if size>0:
                    family = []
                    if len(stack) < size:
                        family = list(stack)
                        for i in family:
                            stack.remove(i)
                    else:
                        for i in range(size):
                            member = rnd.choice(stack)
                            stack.remove(member)
                            family.append(member)
                    list_of_edges[str(edge_ind)]=family
                    edge_ind+=1
        return list_of_edges
    