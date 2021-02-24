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

# The next part is a root-theoretic calculation for computing cut functionals
simples = [var('a'+str(i)) for i in range(1,5)]
roots = [a1, a1+a4, a1+a2+a4, a2, a1+a2+a3+a4,a1+a3+a4, a1+a2+a3+2*a4, a2+a4, a2+a3+a4, a3, a3+a4, a4]
cartanMatrix = matrix([[2,0,0,-1],[0,2,0,-1],[0,0,2,-1],[-1,-1,-1,2]])

def root_to_vec(r):
    return vector([r.coefficient(x) for x in simples])

def reflect(r1,r2):
    """This is the reflection corresponding to the root r1, applied to the root r2."""
    ip = root_to_vec(r1)*cartanMatrix*root_to_vec(r2)
    return (r2 - ip*r1)

simple_objects_dict = {}
for i in range(0,12):
    if i == 0:
        # Start with original list of simples.
        simple_objects_dict[0] = simples
    else:
        # Modify by moving the last element to the front and applying its reflection.
        simple_system = simple_objects_dict[i-1]
        r = simple_system[-1]
        simple_objects_dict[i] = [-r] + [reflect(r,s) for s in simple_system[:-1]]

# Returns the calculation of the functional on the sequence of roots (aka objects).
def get_functional(i):
    if i > len(roots):
        return []
    else:
        eqs = [simple_objects_dict[i][k] == 1 - sgn(k) for k in range(0, len(simple_objects_dict[i]))]
        sol = solve(eqs, simples)
        return [abs(r.substitute(sol[0])) for r in roots]

# Absolute values of the functionals. This is the cut/pair matrix.        
functionals = []
for i in range(0,len(roots)):
    row = get_functional(i)
    functionals = functionals + [row]

for r in roots:
    for s in roots:
        pairing = root_to_vec(r)*cartanMatrix*root_to_vec(s)
        if pairing == 0:
            print (r,s)

# # The following calculation makes a list of all spherical objects in the heart (plus two extra that are not in the heart).
# # This is in order to find maximal subcollections for which any two objects are pairwise parity.
# b1 = [[x,y,z] for x in [s2,t2] for y in [s3,t3] for z in [s4,t4]]
# b2 = [[w,x,y,z] for x in [s2,t2] for y in [s3,t3] for z in [s4,t4] for w in [s1,t1]]
# b3 = [[x,y] for x in [s2,t2] for y in [s3,t3]]
# b4 = [[x,y] for x in [s2,t2] for y in [s4,t4]]
# b5 = [[x,y] for x in [s3,t3] for y in [s4,t4]]
# b6 = [[x] for x in [s2,t2,s3,t3,s4,t4]]
# braids = b1 + b2 + b3 + b4 + b5 + b6

# objects = [(composeAll(x))(P1) for x in braids] + [P1, P2, P3, P4]

# def isParity(A,B):
#     h = hom(A,B).qPolynomial()
#     h1,h2 = h.substitute(q = 1), h.substitute(q = -1)
#     even,odd = h1+h2, h1-h2
#     if even == 0 or odd == 0:
#         return True
#     else: return False

# # Construct a graph with nodes as objects, and edges between two nodes if and only if homs between the two are parity.
# # Then find cliques within this graph.
# import networkx as nx
# parityPairs = [(objects[i],objects[j]) for i in range(0,len(objects)) for j in range(0,i) if isParity(objects[i],objects[j])]
# G = nx.Graph()
# G.add_edges_from(parityPairs)

# # Find maximal cliques
# cliques = list(nx.find_cliques(G))
# M = max([len(x) for x in cliques]) # This ends up being 10.
# maxCliques = [x for x in cliques if len(x) == M]
