from .util import *

import random

# from pysat.solvers import Glucose4
# from pysat.solvers import Glucose3
# from pysat.solvers import Minisat22
from pysat.solvers import MapleChrono
from pysat.examples.musx import MUSX
from pysat.formula import WCNF


# Generating and verifying SAT problems for graphs.
# Returns the clauses, a map from node ix to clause ix, and vice versa.
def getSAT(graph, num_colours):
    G = graph.G # networkx graph
    ps = allPairs(range(num_colours))

    clauses = []
    node_cl_map = {n: [] for n in G.nodes}  # Initilize empty
    cl_node_map = []

    def addCl(c, nodes):
        clauses.append(c)
        cl_node_map.append(nodes)
        for n in nodes:
            node_cl_map[n].append(len(clauses)-1)

    # Mapping from networkx to first colour var per node.
    conv = lambda i: i*num_colours + 1

    # First node is first colour
    # addCl([1], [0])

    # Each node must be exactly one colour
    for n in G.nodes:
        addCl([conv(n) + i for i in range(num_colours)], [n])

        for (a,b) in ps:
            addCl([-(conv(n)+a), -(conv(n)+b)], [n])

    # Edges mustn't be the same colour
    for (a,b) in G.edges:
        for i in range(num_colours):
            addCl([-(conv(a)+i), -(conv(b)+i)], [a,b])

    return clauses, node_cl_map, cl_node_map

# Get constraints for a given list of nodes to not all be the same colour.
# This is used to ensure the central nodes of deGrey's M are not tricolor.
def sameColourConstraint(nodes, num_colours):
    # Mapping from networkx to first colour var per node.
    conv = lambda i: i*num_colours + 1

    # Add constraint on this, no time that
    # 3 of the nodes are the same colour.
    # CNF (n1c1 && n2c1 && n3c1) || (n1c2 && n2c2 && n3c2)
    #  || (n1c3 && n2c3 && n3c3) || (n1c4 && n2c4 && n3c4)
    # == very long...
    # AND of all combinations of the variables from each bracket above ^.

    # Gives lists of (ix, col) to be OR'd
    def rec(col):
        if col < 0:
            return [[]]

        rs = rec(col-1)
        #res = [r + [(n, col)] for n in nodes for r in rs]  # debug
        res = [r + [conv(n) + col] for n in nodes for r in rs]
        return res

    return rec(num_colours-1)

# Say that 2 nodes must have equal colours.
def colourEqualsConstraint(nodeA, nodeB, num_colours):
    conv = lambda i: i*num_colours + 1

    res = []
    a, b = conv(nodeA), conv(nodeB)

    # a should equal b.
    # which is equilivent to (a v !b) ^ (!a v b)
    # cause CNF magic or w/e.
    for c in range(num_colours):
        res.append([(a+c), -(b+c)])
        res.append([-(a+c), (b+c)])

    return res


# Generate minimal graph for unsat colours.
# Optionally include required additional clauses and nodes to keep.
def genMinGraph(graph, num_colours=4, approx=False, required_cl=[], required_nodes=[]):
    cl, n_clM, cl_nM = getSAT(graph, num_colours)

    f = WCNF()
    for c in required_cl:
        f.append(c)
    for c in cl:
        f.append(c, weight=1)
    print("Created formula")

    if approx:
        required_cl = genApproxMinClauses(f)
    else:
        # Calculate MUS
        mu = MUSX(f, verbosity=2)
        required_cl = mu.compute()

    # Map back to graph
    req_nodes = set(required_nodes)
    [req_nodes.update(cl_nM[i-1]) for i in required_cl]
    print("New graph size:", len(req_nodes))
    # Create blacklist and delete
    bk_lst = [n for n in graph.G.nodes if n not in req_nodes]
    [graph.G.remove_node(i) for i in bk_lst]

def optimize(graph, num_colours=4, extract_MUS=False, required_cl=[], required_nodes=[], verbosity=0, shuffle=False):
    """
    Use clever SAT things to reduce the number of nodes in the graph,
    while still maintaining it's UNSAT status.

    extract_MUS runs an additional optimization that will find a minimal
    graph, but it is more expensive.

    required_cl: provide additional clauses that won't be optimized away

    required_nodes: blacklist some nodes from removal. (Typically those in required_cl)

    verbosity: set logging level
    """

    def maybeLog(s, th=1, nl=True):
        if verbosity >= th:
            print(s, end=('\n' if nl==True else ''))

    clauses, n_clM, cl_nM = getSAT(graph, num_colours)


    f = WCNF()
    for c in required_cl:
        f.append(c)
    for c in clauses:
        f.append(c, weight=1)
    maybeLog("Created formula")


    topv = f.nv # Number of vars in basic formula
    # Selectors for each node, each getting a new variable.
    # Note: graph.G.nodes may not be a contigous list after previous removal steps.
    sels = [i+topv for i in graph.G.nodes]
    # vmap = {}   # Map to origional clause

    s = MapleChrono(bootstrap_with=f.hard)

    # Possibly load the graph nodes in in a random order.
    ixs = [i for i in range(0, len(clauses))]
    if shuffle:
        random.shuffle(ixs)
    # For each node, add the relevent clauses
    for i in ixs:
        clause = clauses[i]
        # Nodes involved in this clause
        nodes = cl_nM[i]

        # Assume node is enabled (as negative value).
        # Need both nodes in an edge on to care about rest of clause.
        s.add_clause(clause + [-(n + topv) for n in nodes])

    if not s.solve(assumptions=sels):
        maybeLog("Core extraction\n")
        approx = s.get_core()

    if extract_MUS:
        # Perform additional refinement to get MUC.
        # Attempt to remove each node.
        i = 0
        while i < len(approx):
            # Try ignoring nodes (setting variable to positive),
            # And seeing what happens...
            to_test = approx[:i] + approx[(i + 1):]
            sel, node = approx[i], approx[i] - topv

            maybeLog('c testing node: {0}'.format(node), nl=False)

            if s.solve(assumptions=to_test):
                maybeLog(' -> sat (keeping {0})'.format(node))
                i += 1
            else:
                maybeLog(' -> unsat (removing {0})'.format(node))
                approx = to_test

    # Map back to node ixs, adding to passed required nodes.
    required_nodes = [x - topv for x in approx] + required_nodes

    # Create blacklist and delete from graph
    bk_lst = [n for n in graph.G.nodes if n not in required_nodes]
    print("Removing", len(bk_lst), "nodes.", len(required_nodes), "left.")
    [graph.G.remove_node(i) for i in bk_lst]
    # TODO: Remove the points too...


# Generate approximate minimal clauses for unsat,
# which is returned as the core from the SAT solver.
def genApproxMinClauses(f):
    sels = []  # Selectors
    vmap = {}  # Map to origional clause
    topv = f.nv

    s = MapleChrono(bootstrap_with=f.hard)

    print("Adding clauses")
    # Relaxing soft clauses and adding them to the solver
    for i, c in enumerate(f.soft):
        topv += 1

        sels.append(topv)
        vmap[topv] = i
        s.add_clause(c + [-topv])

    print("Solving")
    if not s.solve(assumptions=sels):
        # get an overapproximation of an MUS
        print("Core extraction (and mapping)")
        approx = [vmap[x] + 1 for x in s.get_core()]
    else:
        print("Wait fuck it's solvable already")
        return

    print("Required clauses:", len(approx))
    return approx

def isColourable(G, num_colours=4, extra_clauses=[]):
    cl, _, _ = getSAT(G, num_colours)
    s = MapleChrono()
    for c in cl + extra_clauses:
        s.add_clause(c)

    return s.solve()
