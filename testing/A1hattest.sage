import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

# * Load files:
load("../complexes.sage")
load("../zigzagalgebra.sage")
load("../zigzagmodules.sage")
load("../braidactions.sage")
load("../HN.sage")

# * The A1-hat quiver and algebra:

a1hatgraph = DiGraph({1: {2: ['a1','a2']}, 2:{1: ['b1','b2']}})
a1hatquiv = make_test(a1hatgraph)
R = a1hatquiv['Z']
R.inject_variables()

F1 = ZigZagModule(R, 1, name="P1")
F2 = ZigZagModule(R, 2, name="P2")

P1 = ProjectiveComplex(R)
P1.addObject(0, F1)

P2 = ProjectiveComplex(R)
P2.addObject(0, F2)

# * Standard twists:

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

def t1(C):
    return t(1, C)

def t2(C):
    return t(2, C)


