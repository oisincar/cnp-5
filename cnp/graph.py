from .util import *
import networkx as nx
import matplotlib.pyplot as plt

import sympy as sp
from sympy import Symbol, I, S
import scipy.spatial as spatial

class Graph():
    def __init__(self, PPV, approx=True):
        """
        Construct a graph from it's points, which should be given
        as a tuple of exact sympy values and approximate python values.
        
        All unit length edges are efficiently generated.

        Set approx=False to use sympy to verify all constructed edges are exactly
        1 unit long. This can be expensive for larger graphs.
        """
        P, PV = PPV
        self.P = P
        self.PV = PV
        self.G = nx.Graph()
        for ix, pV in enumerate(PV):
            self.G.add_node(ix)
            self.G.nodes[ix]['pos'] = (pV.real,pV.imag)

        # Create edges
        psV = csToPnts(PV)
        tree = spatial.cKDTree(psV)
        _accuracy = 0.0000000001
        # dss, nixss = tree.query(psV, max_adjacency, distance_upper_bound=_accuracy+1)
        # nixss = tree.query_ball_point()
        # for ix, (ds, nixs) in enumerate(zip(dss, nixss)):
        #     for i in range(len(ds)): # Each neighbour candidate
        #         if abs(ds[i]-1) < _accuracy:
        #             self.G.add_edge(ix, int(nixs[i]))
                    # TODO: Test for exact measurement.
                   
        pairs = tree.query_pairs(1+_accuracy)
        for (a,b) in pairs:
            if abs(dist(PV[a], PV[b])-1) < _accuracy:
                self.G.add_edge(a, b)
                # Optionally check for exact measurement here.
                
    def checkNetXGraph(self, log_progress=False):
        """
        Check all edges are exactly length 1
        """
        print("Checking:", len(self.G.edges), "edges.")

        for ix, (a,b) in enumerate(self.G.edges):
            c1, c2 = self.P[a], self.P[b]
            if sp.simplify(sp.Abs(c1-c2)-1) != 0:
                print("\n\n:c .. c1:", c1, "\n", c2)
                return False
            if log_progress and ix % (len(self.G.edges)//100) == 0:
                print((ix*100)//len(self.G.edges), "%")

        return True

    def show(self, automate_pos = False):
        """
        Draw the graph.
        Set automate_pos = True to automatically
        determine point locations.
        """
        if automate_pos:
            nx.draw(self.G)
        else:
            nx.draw(self.G, nx.get_node_attributes(self.G, 'pos'), node_size=4)
        plt.show()

def triangle():
    F = nx.Graph()
    F.add_node(0)
    F.add_node(1)
    F.add_node(2)
    F.add_edge(0,1)
    F.add_edge(0,2)
    F.add_edge(1,2)
    return F
