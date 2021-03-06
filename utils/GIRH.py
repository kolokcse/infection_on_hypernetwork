import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from numba import jit
from scipy.spatial.distance import pdist, squareform, cdist
import hypernetx as hnx

def from_incidence_mtx(I_mtx):
    H=hnx.Hypergraph()
    for i in range(I_mtx.shape[1]):
        H.add_edge(hnx.Entity(str(i), elements=[]))
    for i,cols in enumerate(I_mtx.T):
        for v in np.arange(I_mtx.shape[0])[cols==1]:
            H.add_node_to_edge(v, str(i))
    return H 



def _generate_weights(tau=None,
                      expected_weight=None,
                      size=None):
    """Return numpy array of weights"""
    return np.random.random(size)**(-1/(tau-1))*(expected_weight)*(tau-2)/(tau-1)


class GIRH:
    """
    This interpretation of Geometric Inhomogeneous Random Hypergraph based on Geometric Inhomogeneous Random Graph.
    We just made some modification on the code, so we can generate the bipartite graph of the hypergraph as in the
    as a GIRG and generate the hypergraph from the bipartite graph.
    
    Geometric Inhomogeneous Random Graph. 
    Vertices are placed uniformly in a box/torus of length n^(1/d), where n is the number of vertices. 
    Vertices are equipped with an i.i.d. weight following a power law with specified expectation,
    i.e., of the form P(W>w) = cx^{-(\tau-1)}.
    Conditionally on the weights, two vertices x, y are connected with probability 
    p(x,y) = C_1*min(1, C_2*(W_x W_y / |x-y|^d)^alpha) [C2, C1],
    where |x-y| denotes the two-norm in Euclidean space, or on the torus.
    
    Note that there are more general functions p(x,y) described in literature.
    
    Geometric inhomogeneous random graphs
    Karl Bringmann, Ralph Keusch, Johannes Lengler
    Theoretical Computer Science, Volume 760, 14 February 2019, p. 35-54
    https://doi.org/10.1016/j.tcs.2018.08.014
    
    Scale-free Percolation 
    Maria Deijfen, Remco van der Hofstad, Gerard Hooghiemstra
    Annales de l'I.H.P. ProbabilitÃƒÆ’Ã‚Â©s et statistiques, Volume 49 (2013) no. 3, p. 817-838 
    https://doi.org/10.1214/12-AIHP480
    """
    def __init__(self,
                 vertex_size,
                 number_of_hyperedges,
                 alpha,
                 tau,
                 expected_node_weight,
                 expected_edge_weight,
                 C_1=1,
                 C_2=1,
                 on_torus=True,
                 dimension=2):
        self.vertex_size = vertex_size
        self.number_of_hyperedges = number_of_hyperedges
        self.alpha = alpha
        self.tau = tau
        self.expected_node_weight = expected_node_weight
        self.expected_edge_weight = expected_edge_weight
        self.dimension = dimension
        self.C_1 = C_1
        self.C_2 = C_2
        self.vertex_weights = _generate_weights(tau=tau,
                                        expected_weight=expected_node_weight,
                                        size=vertex_size)
        self.edge_weights = _generate_weights(tau=tau,
                                        expected_weight=expected_edge_weight,
                                        size=number_of_hyperedges)
        self.vertex_locations = self.generate_locations(vertex_size)
        self.edge_locations = self.generate_locations(number_of_hyperedges)
        
        self.hypergraph = self.generate_hypergraph(vertex_weights=self.vertex_weights,
                                         edge_weights=self.edge_weights,
                                         vertex_locations=self.vertex_locations,
                                         edge_locations=self.edge_locations,
                                         on_torus=on_torus)

    def _compute_isolated_star_attrs(self):
        isolated_star_dict = dict([(node, _isolated_star_number(self.graph, node))
                                   for node in self.graph.nodes()])
        nx.set_node_attributes(self.graph, isolated_star_dict, name="isolated_star")


    @property
    def torus_width(self):
        return self._torus_width(self.vertex_size, self.dimension)

    def _torus_width(self, vertex_size, dimension):
        return vertex_size**(1/dimension)

    def generate_locations(self,
                           size=None,
                           dimension=None):
        """Returns (size*dimension) numpy array with points on the d-dimensional torus
        of length size^(1/dimension)"""
        size = size or self.vertex_size
        dimension = dimension or self.dimension
        return (np.random.random(size=(size, dimension))
                * self._torus_width(size, dimension))

    def torus_distance(self,
                       point1,
                       point2,
                       torus_width=None):
        torus_width = torus_width or self.torus_width
        abs_diff = np.abs(point1 - point2)
        torus_diff = np.minimum(torus_width - abs_diff, abs_diff)
        return np.sum(torus_diff**2)**0.5
    
    def generate_incidence_matrix(self,
                                  vertex_weights,
                                  edge_weights,
                                  vertex_locations,
                                  edge_locations,
                                  alpha=None,
                                  torus_width=None,
                                  on_torus=True):
        alpha = alpha or self.alpha
        torus_width = torus_width or self.torus_width
        vertex_size = len(vertex_weights)
        number_of_edges = len(edge_weights)
        #weightsum = np.sum(vertex_weights)+np.sum(edge_weights)

        dimension = vertex_locations.shape[1]
        weights_multiplied = np.dot(vertex_weights.reshape(len(vertex_weights), 1),
                                    edge_weights.reshape(1, len(edge_weights)))
        
        distance_metric = self.torus_distance if on_torus else "euclidean"
        
        distances_pairwise = cdist(vertex_locations, edge_locations,
                                   distance_metric)
        
        #np.fill_diagonal(distances_pairwise, val=1)
        #p(x,y) = C_1*min(1, C_2*(W_x W_y / |x-y|^d)^alpha) [C2, C1],
        #edge_probabilities = 1 - np.exp(-(weights_multiplied / (distances_pairwise**dimension))**alpha)
        edge_probabilities = self.C_1*np.minimum(
            np.ones(weights_multiplied.shape),
            self.C_2*(weights_multiplied / ((distances_pairwise**dimension)))**alpha)
            
            
        edge_rvs_uniform = np.random.random(size=(vertex_size, number_of_edges))
        incidence_matrix = (edge_rvs_uniform < edge_probabilities).astype(np.int)
        return incidence_matrix

    def generate_hypergraph(self,
                       vertex_weights,
                       edge_weights,
                       vertex_locations,
                       edge_locations,
                       alpha=None,
                       torus_width=None,
                       on_torus=True):
        vertex_weight_dict = dict(zip(range(len(vertex_weights)), vertex_weights))
        vertex_position_dict = dict(zip(range(len(vertex_locations)), vertex_locations))
        edge_weight_dict = dict(zip(range(len(edge_weights)), edge_weights))
        edge_position_dict = dict(zip(range(len(edge_locations)), edge_locations))
        I_mtx = self.generate_incidence_matrix(vertex_weights=vertex_weights,
                                                edge_weights=edge_weights,
                                                vertex_locations=vertex_locations,
                                                edge_locations=edge_locations,
                                                alpha=alpha,
                                                torus_width=torus_width,
                                                on_torus=on_torus)
        H=from_incidence_mtx(I_mtx)
        #nx.set_node_attributes(graph, weight_dict, name="weight")
        #nx.set_node_attributes(graph, position_dict, name="position")
        return H


    def draw(self, pos=None):
        pos={}
        for node_i,node in enumerate(self.hypergraph.nodes):
            pos[node_i]=tuple(self.vertex_locations[node_i])
        for edge_i,edge in enumerate(self.hypergraph.edges):
            pos[str(edge_i)]=tuple(self.edge_locations[edge_i])
        hnx.draw(H, pos=pos)