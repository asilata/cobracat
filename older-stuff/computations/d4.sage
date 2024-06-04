load("../complexes.sage")
load("../zigzagalgebra.sage")
load("../zigzagmodules.sage")
load("../braidactions.sage")

d4 = make_test(d4graph)
R = d4['Z']
R.inject_variables()
F = ZigZagModule(R, 1, name = "P1")
G = ZigZagModule(R, 2, name = "P2")
H = ZigZagModule(R, 3, name = "P3")
J = ZigZagModule(R, 4, name = "P4")

P1 = ProjectiveComplex(R)
P1.addObject(0, F)

P2 = ProjectiveComplex(R)
P2.addObject(0, G)

P3 = ProjectiveComplex(R)
P3.addObject(0, H)

P4 = ProjectiveComplex(R)
P4.addObject(0, J)


def s(i, C):
    D = sigma(R, i, C)
    D.minimize()
    return D

def t(i, C):
    D = sigmaInverse(R, i, C)
    D.minimize()
    return D

def s1(C):
    return s(1, C)

def s2(C):
    return s(2, C)

def s3(C):
    return s(3, C)

def s4(C):
    return s(4, C)

def t1(C):
    return t(1, C)

def t2(C):
    return t(2, C)

def t3(C):
    return t(3, C)

def t4(C):
    return t(4, C)


heart = [P1, P2, P3, P4]

# The following calculation makes a list of all spherical objects in the heart (plus two extra that are not in the heart).
# This is in order to find maximal subcollections for which any two objects are pairwise parity.
b1 = [[x,y,z] for x in [s2,t2] for y in [s3,t3] for z in [s4,t4]]
b2 = [[w,x,y,z] for x in [s2,t2] for y in [s3,t3] for z in [s4,t4] for w in [s1,t1]]
b3 = [[x,y] for x in [s2,t2] for y in [s3,t3]]
b4 = [[x,y] for x in [s2,t2] for y in [s4,t4]]
b5 = [[x,y] for x in [s3,t3] for y in [s4,t4]]
b6 = [[x] for x in [s2,t2,s3,t3,s4,t4]]
braids = b1 + b2 + b3 + b4 + b5 + b6

objects = [(composeAll(x))(P1) for x in braids] + [P1, P2, P3, P4]

def isParity(A,B):
    h = hom(A,B).qPolynomial()
    h1,h2 = h.substitute(q = 1), h.substitute(q = -1)
    even,odd = h1+h2, h1-h2
    if even == 0 or odd == 0:
        return True
    else: return False

# Construct a graph with nodes as objects, and edges between two nodes if and only if homs between the two are parity.
# Then find cliques within this graph.
import networkx as nx
parityPairs = [(objects[i],objects[j]) for i in range(0,len(objects)) for j in range(0,i) if isParity(objects[i],objects[j])]
G = nx.Graph()
G.add_edges_from(parityPairs)

# Find maximal cliques
cliques = list(nx.find_cliques(G))
M = max([len(x) for x in cliques]) # This ends up being 10.
maxCliques = [x for x in cliques if len(x) == M]
