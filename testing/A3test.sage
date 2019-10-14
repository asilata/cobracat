import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

# * Load files:
load("../complexes.sage")
load("../zigzagalgebra.sage")
load("../zigzagmodules.sage")
load("../braidactions.sage")
load("../HN.sage")

# * The A3 quiver and algebra:

a3graph = DiGraph({1: {2: 'a1'}, 2:{1: 'a2', 3:'b1'}, 3:{2:'b2'}})
a3quiv = make_test(a3graph)
R = a3quiv['Z']
R.inject_variables()

# * Standard objects:
F1 = ZigZagModule(R, 1, name="P1")
F2 = ZigZagModule(R, 2, name="P2")
F3 = ZigZagModule(R, 3, name="P3")

P1 = ProjectiveComplex(R)
P1.addObject(0, F1)

P2 = ProjectiveComplex(R)
P2.addObject(0, F2)

P3 = ProjectiveComplex(R)
P3.addObject(0, F3)

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

def s3(C):
    return s(3, C)

def t1(C):
    return t(1, C)

def t2(C):
    return t(2, C)

def t3(C):
    return t(3, C)

# * HN Filtrations (convex square stability condition):
# Let us go through the charts

HN(P2, stab)
HN(s3(P2), stab)
HN(sz(s3(P2)), stab)
HN(s1(sz(s3(P2))), stab)
HN(s2(s1(sz(s3(P2)))), stab)    
HN(t3(s2(s1(sz(s3(P2))))),stab)
HN(s2(t3(s2(s1(sz(s3(P2)))))),stab)

# *

# ** Stable objects:
x = s1(P2)
y = s2(P3)
z = s1(y)
# ** Stable twists:
sx = composeAll([s1,s2,t1])
sy = composeAll([s2,s3,t2])
sx = composeAll([s1,s2,t1])
sz = composeAll([s1,s2,s3,t2,t1])
# ** Stability condition:
stab = [P3, s2(P3), s1(s2(P3)), P2, s1(P2),P1]

# * HN Filtrations (concave square stability condition):
# ** Stable objects:
X = s2(P1)
Y = s2(P3)
Z = s1(s2(P3))
sx = composeAll([s2,s1,t2])
sy = composeAll([s2,s3,t2])
sz = composeAll([s1,s2,s3,t2,t1])

# ** Stability condition:
stab = [P3, Y, Z, P1, X, P2]

# ** Charts:
# Guess by applying standard twists to all stable objects
[HN(s1(o),stab) for o in stab]
# [[[0]: P3<0> :[0]],
#  [[-2]: P1<-2> → P2<-1> → P3<0> :[0]],
#  [[-2]: P1<-2> → P2<-1> → P3<0> :[0], [-3]: P1<-4> :[-3]],
#  [[-1]: P1<-2> :[-1]],
#  [[-1]: P2<-1> :[-1]],
#  [[0]: P2<0> :[0], [-1]: P1<-1> :[-1]]]
[HN(t1(o),stab) for o in stab]
# [[[0]: P3<0> :[0]],
#  [[0]: P1<0> :[0], [-1]: P2<-1> → P3<0> :[0]],
#  [[-1]: P2<-1> → P3<0> :[0]],
#  [[1]: P1<2> :[1]],
#  [[1]: P1<2> :[1], [-1]: P2<-1> → P1<0> :[0]],
#  [[0]: P2<0> → P1<1> :[1]]]
[HN(s2(o),stab) for o in stab]
# [[[-1]: P2<-1> → P3<0> :[0]],
#  [[-1]: P2<-1> → P3<0> :[0], [-2]: P2<-3> :[-2]],
#  [[-2]: P1<-2> → P2<-1> → P3<0> :[0]],
#  [[-1]: P2<-1> → P1<0> :[0]],
#  [[-1]: P2<-1> → P1<0> :[0], [-2]: P2<-3> :[-2]],
#  [[-1]: P2<-2> :[-1]]]
[HN(t2(o),stab) for o in stab]
# [[[1]: P2<1> :[1], [0]: P3<0> :[0]],
#  [[0]: P3<0> :[0]],
#  [[-2]: P1<-2> → P2<-1> → P3<0> :[0]],
#  [[1]: P2<1> :[1], [0]: P1<0> :[0]],
#  [[0]: P1<0> :[0]],
#  [[1]: P2<2> :[1]]]
[HN(s3(o),stab) for o in stab]
# [[[-1]: P3<-2> :[-1]],
#  [[-1]: P2<-1> :[-1]],
#  [[-1]: P2<-1> :[-1], [-2]: P1<-2> :[-2]],
#  [[0]: P1<0> :[0]],
#  [[-1]: P2<-1> → P1<0> :[0], [-2]: P3<-2> :[-2]],
#  [[0]: P2<0> :[0], [-1]: P3<-1> :[-1]]]

# These seem to give all the charts


