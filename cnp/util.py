import math
import cmath

def dist(c1, c2):
    """
    Distance in the complex plane
    """
    c = c2 - c1
    return math.sqrt(c.imag*c.imag + c.real*c.real)

def allPairs(lst):
    """
    All pairs of a list.
    """
    if len(lst) <= 1:
        return []
    
    h, *t = lst
    return [(h,x) for x in t] + allPairs(t)

def allSumPairs(lst):
    """
    Sums of all pairs in a list.
    """
    return [a+b for (a,b) in allPairs(lst)]

def allProd(ws):
    """
    Products of all subsets of a list.
    """
    if len(ws) == 0:
        return []
    
    h, *t = ws
    rec = allProd(t)
    return [h] + rec + [h * k for k in rec]

def cont(pp, lst, approx=True):
    """
    Check if a point is contained in a list (naively)
    Set approx=False to check sympy values simplify to exactly 0.
    """
    a,aV = pp

    for (b, bV) in lst:
        if aV == bV:
            return True

        if dist(aV, bV) < 0.000000001:
            # To speed up code... Disable exact checking.
            # Check formulas simplify to exactly 0
            if approx or sp.simplify(sp.Abs(a-b)) == 0:
                return True

    return False

def csToPnts(cs):
    """
    Convert complex to coord
    """
    return [(c.real, c.imag) for c in cs]
