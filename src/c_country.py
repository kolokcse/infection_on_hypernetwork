from subprocess import Popen, STDOUT, PIPE

import csv
import json
import os
import sys
import copy
import numpy as np
import networkx as nx
import hypernetx as hnx
import multiprocessing
from utils.hypergraph_generator import get_hypergraph
import matplotlib.pyplot as plt

global shm

def iter_by2(iterable):
    return zip(*[iter(iterable)]*2)

def get_Chi(aggregated, cities):
    ratioInf=aggregated[:,2]/(np.sum(aggregated, axis = 1))
    
    xi = cities[:,:,2]
    mi = np.sum(cities,axis=2)
    
    sumChi2 = np.sum((xi-mi)**2/mi, axis=1)
    
    return sumChi2

class C_Country:
    def __init__(self, hypergraph):
        # Parameters:
        #     graph      : network of cities
        self.hypergraph = hypergraph
        self.node_num = self.hypergraph.number_of_nodes()
        self.edge_num = self.hypergraph.number_of_edges()
        self.contactNumber = None

    def get_str_hypergraph(self):
        if self.hypergraph.isstatic:
            edge_neighbors = lambda node: list(self.hypergraph.edges.memberships[node].keys())
        else:
            edge_neighbors = lambda node: list(self.hypergraph.nodes[node].memberships.keys())

        # === City data ===
        str_graph = "{} {}\n".format(self.node_num, self.edge_num)
        for node in self.hypergraph.nodes:
            neighs = [str(k) for k in edge_neighbors(node)]
            line = "{} {}\n".format(len(neighs), " ".join(neighs))
            str_graph += line
        
        for edge in range(self.edge_num):
            neighs = [str(k) for k in list(self.hypergraph.incidence_dict[str(edge)])]
            line = "{} {}\n".format(len(neighs), " ".join(neighs))
            str_graph += line
        return str_graph

    
    @staticmethod
    def set_shm(job_num):
        job_count = multiprocessing.Value("i", 0)
        global shm
        shm = (job_count, job_sum)
    
    @staticmethod
    def run(args, str_hypergraph, node_num, job_count, lock, inf_nodes, verbose=False, log=False):
        # === Infect ===
        # === Inf cities ===
        str_inf_nodes = "{}\n".format(len(inf_nodes))
        for inf_node in inf_nodes:
            str_inf_nodes += "{}\n".format(inf_node)
        
        # === Agent data ===
        # TODO: agents should be initialized in C++ 
        str_args = [str(item) for pair in args.items() for item in pair]
        p = Popen([os.path.dirname(os.path.abspath(__file__))+'/bin/main']+ str_args,
                  stdout=PIPE, stdin=PIPE, stderr=STDOUT, bufsize=1, universal_newlines=True)

        out,err = p.communicate(str_hypergraph + str_inf_nodes)

        history = []
        edges = []
        
        if(verbose):
            for line,line2 in iter_by2(out.split('\n')[:-1]):
                history.append([int(a) for a in line[:-1].split(" ")])
                edges.append([[int(s) for s in edge[:-1].split(" ")] for edge in line2[:-1].split(";")])
        else:
            for line in out.split('\n')[:-1]:
                history.append([int(a) for a in line[:-1].split(" ")])        
                
        history = np.array(history[:])
        edges = np.array(edges)
        
        
        with lock:
            job_count[0]+=1
            print('\r {}/{}'.format(job_count[0], job_count[1]), end='', flush=True)
            
        return history,edges
   