# Set of functions for generating the coordinates of
# interesting graphs, including those from the origional paper.

# All functions return pairs of points. Those in exact sympy form,
# and those in approximate python notation.

import math
import cmath

import sympy as sp
from sympy import Symbol, I, S
from sympy.parsing.sympy_parser import parse_expr

import scipy.spatial as spatial
from .util import *

# Error for floating point calculations.
_accuracy = 0.00000001

# As first mentioned here:
# https://dustingmixon.wordpress.com/2018/05/01/polymath16-third-thread-is-6-chromatic-within-reach/
# Minimal graphs discovered so far are subsets of Z[w1, ..]
def w(tV):
    return cmath.exp(math.acos(1 - 1/(2*tV)) * 1j)
def w_(t):
    return sp.exp(I*sp.acos(1 - 1/(2*S(t))))

# Shift graph G by offset c
def shifted(PPV, c):
    P,PV = PPV
    cV = complex(c)
    return ([p+c for p in P], # Possibly simplify here?
            [pV+cV for pV in PV])

# Rotate about point cOff (complex), by angle ang (rad)
def rotated(GGV, offset, angle, reEval=False):
    G, GV = GGV
    cAng = sp.simplify(sp.cos(angle) + I*sp.sin(angle))
    #G2 = [sp.simplify(cOff + (g-cOff)*cAng) for g in G]
    G2 = [offset + (g-offset)*cAng for g in G]

    if reEval:
        GV2 = [complex(g) for g in G2]
    else:
        cOV,cAV = complex(offset), complex(cAng)
        GV2 = [cOV + (gv-cOV)*cAV for gv in GV]

    return (G2, GV2)

def from_complex(ps):
    return ([p for p in ps],
            [complex(p) for p in ps])

def empty():
    return ([],[])

# Add sympy point to some other points.
def addP(pSet, p):
    # TODO: Check if it's already here...

    pSet[0].append(p)
    pSet[1].append(complex(p))
    return pSet


# Add graph H to G. G should be larger.
def add(GGV, HHV):
    G, GV = GGV
    H, HV = HHV

    if G == []:
        return (H, HV)
    if H == []:
        return (G, GV)

    # Find list of very nearby points to each new point
    tree = spatial.cKDTree(csToPnts(GV))
    # bList = tree.query_ball_point(csToPnts(HV), _accuracy)

    tree2 = spatial.cKDTree(csToPnts(HV))
    bList = tree2.query_ball_tree(tree, _accuracy)

    for ix, (h, hV) in enumerate(zip(H, HV)):
        if len(bList[ix]) == 0:
            G.append(h)
            GV.append(hV)

#         else:  # Check exact math (Untested!!)
#             found = False
#             for nIx in bList[ix]:
#                 if sp.simplify(sp.Abs(a-b)) == 0:
#                     found = True
#                     break
#             if found = False:
#                 print("Ah fuck, didn't simplify")

    return (G, GV)

# Hexagonal grid piece w/ origin at ix 0.
def makeH():
    w = sp.simplify((1 + I * sp.sqrt(3))/2)
    H = [0*I] + [sp.simplify(w**i) for i in range(6)]
    return (H, [complex(h) for h in H])

# Grid of hexagons (or any grid passed to it) (Pass M here)
def makeJ(core=None):
    # Grab offsets from H
    H, _ = makeH()

    if core == None:
        MMV = makeMGraph()
    else:
        MMV = core

    # Find unique offset points
    centers = []
    [centers.append(pp)
        for pp in [(p, complex(p)) for p in [0*I] + allSumPairs(H)]
        if not cont(pp, centers)]
    print("Creating:", len(centers))
    JJV = ([],[])
    # Sum pairs of offsets to get the coords of the outer hexagons.
    for (c,_) in centers:
        MMV_sh = shifted(MMV, c)
        JJV = add(JJV, MMV_sh)

    return JJV

# Two Js, rotated.
def makeK(core=None):
    J = makeJ(core)
    print("J made")
    Jr = rotated(J, offset=0*I, angle=2*sp.asin(S(1)/4))
    print("Adding")
    return add(J, Jr)

# Two Ks, rotated.
def makeL(core=None):
    K = makeK(core)
    print("K made")
    Kr = rotated(K, offset=-2 + 0*I, angle=2*sp.asin(S(1)/8))
    print("Adding")
    return add(K, Kr)

# DeGrey's origional N with 52 copys of M.
def makeN():
    return makeL(makeM())

# TODO clean with add (, rotated?)
def makeM():
    # Generate using angles
    # a, b = sp.symbols('a b')
    # expr =
    # angs = [a * sp.asin(sp.sqrt(3)/2) + b * sp.asin(1/sp.sqrt(12))
    #         for b in range(-2, 3) for a in range(0,6) ]
    # ps = [sp.simplify(sp.cos(a) + I*sp.sin(a)) for a in angs]

    # allPs = []
    # [allPs.append(p)
    #      for p in allSumPairs(ps)
    #      if not cont((p, complex(p)), allPs)
    #         and abs(complex(p)) <= math.sqrt(3) + _accuracy]

    # Vectorspace of graph is H rotated different amounts.
    H = makeH()
    Hs = ([],[])
    angs = [a * sp.asin(1/sp.sqrt(12)) for a in range(-2, 3)]
    for a in angs:
        Hs = add(Hs, rotated(H, offset=0, angle=a))

    # Construct all pairs in the vector space.
    # Starting with H means this will make up verts 0-6 in the final result.
    Ht = ([],[])
    print(len(Hs[0]))
    for off in Hs[0]:
        Ht = add(Ht, shifted(Hs, off))

    # Filter ones that are too long, and convert back to tuple.
    hs = [h for h in Ht[0] if dist(complex(h),0) <= math.sqrt(3) + _accuracy]
    Ht = (hs, [complex(h) for h in hs])

    # Finally, construct one of these at each vertice in H.
    M = makeH()
    for off in makeH()[0]:
        M = add(M, shifted(Ht, off))

    return M

    # Now with all vectors, create all combinations sums of 1 or 2

    # Adding origin and ps first makes central H graph show up first in result.
    # vectPs = [(v, complex(v)) for v in [0*I] + ps + allSumPairs(ps)]

    # Construct central graph, making sure the origin is first.
    # Verts 0-6 are the central H.
    # allPs = []
    # [allPs.append(pp)
    #      for pp in vectPs
    #      if not cont(pp, allPs)
    #         and dist(pp[1], 0) <= math.sqrt(3) + _accuracy]

    # # Now!!
    # # Construct a new graph M with 7 duplicates of the previous
    # w = (1 + I * sp.sqrt(3))/2
    # offsets = [sp.simplify(w**i) for i in range(6)]

    # M = [p for (p,_) in allPs]    # Exact points
    # MV = [pV for (_,pV) in allPs] # Approx points

    # for off in offsets:
    #     offV = complex(off)

    #     # Find list of very nearby points to each new point
    #     tree = spatial.cKDTree(csToPnts(MV))
    #     q = [pV + offV for (_, pV) in allPs] # Coords of new copy
    #     bList = tree.query_ball_point(csToPnts(q), _accuracy)

    #     #print(len(allPs), len(bSet))
    #     for ix, (p, pV) in enumerate(allPs):
    #         if len(bList[ix]) == 0:
    #             M.append(p + off)
    #             MV.append(pV + offV)
    #         else:  # Check exact math (Untested!!)
    #             found = False
    #             for nIx in bList[ix]:
    #                 if sp.simplify(sp.Abs(a-b)) == 0:
    #                     found = True
    #                     break
    #             if found = False:
    #                 print("Ah fuck, didn't simplify")


    print("Full graph with", len(M), "vertices")
    return (M, MV)


def fromFile(f_name):
    with open(f_name) as f:
        points_str = f.read().splitlines()

    num_lines = len(points_str)
    i = 0

    def readExpr(s):
        s = s.replace('Sqrt', 'sqrt').replace('[', '(').replace(']', ')')
        return parse_expr(s, evaluate=True)

    def readPoint(ix, pstr):
        # Print progress
        if num_lines > 1000 and ((ix-1)*100) // num_lines != (ix*100) // num_lines:
            print((ix*100)//num_lines, '%')

        # Drop first and last char {/}
        pstr = pstr[1:-1]
        cs = pstr.split(',')
        return readExpr(cs[0]) + readExpr(cs[1]) * I

    return [readPoint(ix, pstr) for ix, pstr in enumerate(points_str)]
