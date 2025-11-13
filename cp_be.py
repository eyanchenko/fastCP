import cpnet
import numpy as np
import networkx as nx

# Code to detect CP labels using cpnet package
# Input: Adjacency matrix, A
# Output: CP labels, C

def cp_be(A, rns=10):
    
    algorithm = cpnet.BE(num_runs=rns)
    G = nx.from_numpy_array(A)
    algorithm.detect(G)
    x = algorithm.get_coreness()
    
    return list(x.values())
