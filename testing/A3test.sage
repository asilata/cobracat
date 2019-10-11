# * Load files
load("../complexes.sage")
load("../zigzagalgebra.sage")
load("../zigzagmodules.sage")
load("../braidactions.sage")
load("../HN.sage")

# * The A3 quiver and algebra

a3graph = DiGraph({1: {2: 'a1'}, 2:{1: 'a2', 3:'b1'}, 3:{2:'b2'}})
a3quiv = make_test(a3graph)
R = a3quiv['Z']
R.inject_variables()

# * Standard objects
F1 = ZigZagModule(R, 1, name="P1")
F2 = ZigZagModule(R, 2, name="P2")
F3 = ZigZagModule(R, 3, name="P3")

P1 = ProjectiveComplex(R)
P1.addObject(0, F1)

P2 = ProjectiveComplex(R)
P2.addObject(0, F2)

P3 = ProjectiveComplex(R)
P3.addObject(0, F3)

# * Standard twists

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

def t1(C):
    return t(1, C)

def t2(C):
    return t(2, C)

def t3(C):
    return t(3, C)

# * Tests
x = s1(P2)
y = s2(P3)
z = s1(y)

sx = composeAll([s1,s2,t1])
sy = composeAll([s2,s3,t2])
sx = composeAll([s1,s2,t1])
sz = composeAll([s1,s2,s3,t2,t1])
stab = [P3, s2(P3), s1(s2(P3)), P2, s1(P2),P1]

# * HN Filtrations
# Let us go through the charts

HN(P2, stab)
HN(s3(P2), stab)
HN(sz(s3(P2)), stab)
HN(s1(sz(s3(P2))), stab)
HN(s2(s1(sz(s3(P2)))), stab)    
HN(t3(s2(s1(sz(s3(P2))))),stab)
HN(s2(t3(s2(s1(sz(s3(P2)))))),stab)
