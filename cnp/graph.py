from .util import *
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

import sympy as sp
from sympy import Symbol, I, S
import scipy.spatial as spatial

class Graph():
    def __init__(self, PPV, approx=True, large=False):
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
        _accuracy = 0.0000000001
        # dss, nixss = tree.query(psV, max_adjacency, distance_upper_bound=_accuracy+1)
        # nixss = tree.query_ball_point()
        # for ix, (ds, nixs) in enumerate(zip(dss, nixss)):
        #     for i in range(len(ds)): # Each neighbour candidate
        #         if abs(ds[i]-1) < _accuracy:
        #             self.G.add_edge(ix, int(nixs[i]))
                    # TODO: Test for exact measurement.
        psV = csToPnts(PV)

        if not large:
            tree = spatial.cKDTree(psV)

            pairs = tree.query_pairs(1+_accuracy)
            for (a,b) in pairs:
                if abs(dist(PV[a], PV[b])-1) < _accuracy:
                    self.G.add_edge(a, b)
                    # Optionally check for exact measurement here.

        else:
            seg_size = 100
            # Enumerate segments. Start at 1 so we've already added a single point to the graph.
            # Then add the rest in increments.
            for start_ix in range(2, len(PV), seg_size):
                print("building...", (start_ix*100) // len(PV), "%", end='\r')
                end_ix = min(start_ix + seg_size, len(PV))

                # Trees for points up to now, and points in the segment we're considering.
                P_seg = psV[start_ix:end_ix]
                T_seg = spatial.cKDTree(P_seg)

                # find all connections within these points.
                # TODO: Just take difference between pairs at 1+accuracy and 1-accuracy
                pairs1 = T_seg.query_pairs(1+_accuracy)
                pairs2 = T_seg.query_pairs(1-_accuracy)
                ps = list(set(pairs1) - set(pairs2))
                for (a,b) in ps:
                    # if abs(dist(PV[a+start_ix], PV[b+start_ix])-1) < _accuracy:
                    self.G.add_edge(a+start_ix, b+start_ix)

                # Find all connections between these points and the ones we're already seen.
                P_d = psV[:start_ix]
                T_d = spatial.cKDTree(P_d)

                def ptsForQuery(res):
                    a = []
                    for ix, lst in enumerate(res):
                        a.extend([(ix,l) for l in lst])

                    return a

                pairs1 = ptsForQuery(T_d.query_ball_tree(T_seg, 1+_accuracy))
                pairs2 = ptsForQuery(T_d.query_ball_tree(T_seg, 1-_accuracy))
                # print(pairs1, pairs2)
                ps = list(set(pairs1) - set(pairs2))

                for (a,b) in ps:
                    # a is the index of the point we've already seen.
                    self.G.add_edge(a, b + start_ix)

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

    def show(self, automate_pos=False, with_labels=False):
        """
        Draw the graph.
        Set automate_pos = True to automatically
        determine point locations.
        """
        if automate_pos:
            nx.draw(self.G, with_labels=with_labels)
        else:
            nx.draw(self.G, nx.get_node_attributes(self.G, 'pos'), node_size=4, with_labels=with_labels, font_size=25)
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
