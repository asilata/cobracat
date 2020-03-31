
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

# HN(P2, stab)
# HN(s3(P2), stab)
# HN(sz(s3(P2)), stab)
# HN(s1(sz(s3(P2))), stab)
# HN(s2(s1(sz(s3(P2)))), stab)    
# HN(t3(s2(s1(sz(s3(P2))))),stab)
# HN(s2(t3(s2(s1(sz(s3(P2)))))),stab)
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
# ** Unimodality
def HNMultiplicities(ob, stab):
    hnfiltration = HN(ob,stab)
    mults = {s:0 for s in stab}
    for h in hnfiltration:
        for s in stab:
            qPoly = hom(h,s).qPolynomial()
            if qPoly.substitute(q=1) == 2 and qPoly.substitute(q=-1).abs() == 2:
                mults[s] = mults[s]+1
    return mults

## ** Local minimum
beta = composeAll([s1,sx,sy])
def HNMultiplicitiesP(beta, stab):
    mults = {s: 0 for s in stab}
    for s in stab:
        hns = HNMultiplicities(beta(s), stab)
        for m in mults.keys():
            mults[m] = mults[m] + hns[m]
    return mults

# * HN Filtrations (concave square stability condition):
# ** Stable objects:
X = s2(P1)
Y = s2(P3)
Z = s1(s2(P3))
sx = composeAll([s2,s1,t2])
tx = composeAll([s2,t1,t2])
sy = composeAll([s2,s3,t2])
ty = composeAll([s2,t3,t2])
sz = composeAll([s1,s2,s3,t2,t1])
tz = composeAll([s1,s2,t3,t2,t1])

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
[HN(t3(o),stab) for o in stab]
# [[[1]: P3<2> :[1]],
#  [[1]: P3<2> :[1], [-1]: P2<-1> → P3<0> :[0]],
#  [[1]: P3<2> :[1], [-2]: P1<-2> → P2<-1> → P3<0> :[0]],
#  [[0]: P1<0> :[0]],
#  [[0]: P1<0> :[0], [-1]: P2<-1> → P3<0> :[0]],
#  [[0]: P2<0> → P3<1> :[1]]]
[HN(sx(o),stab) for o in stab]
# [[[0]: P1<0> :[0], [-1]: P2<-1> → P3<0> :[0]],
#  [[-1]: P2<-1> → P3<0> :[0]],
#  [[-2]: P1<-2> → P2<-1> → P3<0> :[0], [-4]: P2<-5> → P1<-4> :[-3]],
#  [[0]: P1<0> :[0], [-2]: P2<-3> → P1<-2> :[-1]],
#  [[-2]: P2<-3> → P1<-2> :[-1]],
#  [[0]: P1<1> :[0]]]
[HN(tx(o),stab) for o in stab]
# [[[1]: P2<1> → P1<2> :[2], [0]: P3<0> :[0]],
#  [[-1]: P2<-1> → P3<0> :[0]],
#  [[-1]: P2<-1> → P3<0> :[0], [-2]: P2<-3> :[-2]],
#  [[0]: P2<-1> :[0]],
#  [[0]: P2<1> → P1<2> :[1]],
#  [[1]: P2<2> → P1<3> :[2], [0]: P2<0> :[0]]]
[HN(sy(o),stab) for o in stab]
# [[[0]: P3<0> :[0], [-2]: P2<-3> → P3<-2> :[-1]],
#  [[-2]: P2<-3> → P3<-2> :[-1]],
#  [[-2]: P1<-2> :[-2]],
#  [[0]: P1<0> :[0], [-1]: P2<-1> → P3<0> :[0]],
#  [[-1]: P2<-1> → P1<0> :[0]],
#  [[0]: P3<1> :[0]]]
[HN(ty(o),stab) for o in stab]
# [[[0]: P2<-1> :[0]],
#  [[0]: P2<1> → P3<2> :[1]],
#  [[0]: P2<1> → P3<2> :[1], [-2]: P1<-2> → P2<-1> → P3<0> :[0]],
#  [[0]: P1<0> → P2<1> → P3<2> :[2]],
#  [[-1]: P2<-1> → P1<0> :[0]],
#  [[1]: P2<2> → P3<3> :[2], [0]: P2<0> :[0]]]
[HN(sz(o),stab) for o in stab]
# [[[0]: P3<0> :[0], [-3]: P1<-4> → P2<-3> → P3<-2> :[-1]],
#  [[-1]: P2<-1> → P3<0> :[0], [-3]: P1<-4> → P2<-3> → P3<-2> :[-1]],
#  [[-3]: P1<-4> → P2<-3> → P3<-2> :[-1]],
#  [[0]: P2<1> → P3<2> :[1]],
#  [[0]: P2<1> → P3<2> :[1], [-1]: P2<-1> :[-1]],
#  [[0]: P2<0> :[0]]]
[HN(tz(o),stab) for o in stab]
# [[[0]: P2<-1> :[0], [-1]: P1<-2> :[-1]],
#  [[-1]: P1<-2> :[-1]],
#  [[-1]: P1<0> → P2<1> → P3<2> :[1]],
#  [[1]: P1<2> → P2<3> → P3<4> :[3], [0]: P1<0> :[0]],
#  [[1]: P1<2> → P2<3> → P3<4> :[3], [-1]: P2<-1> → P1<0> :[0]],
#  [[0]: P2<0> :[0]]]

# These seem to give all the charts
# ** Transitions:
# *** Good maps for P1 -> P2 extension:
[[HN(twist(P1), stab), HN(twist(P2), stab), HN(twist(s1(P2)), stab)] for twist in [s1,t1,s2,t2,s3,t3,sx,tx, sy, ty, sz, tz]]
# Good: s1, t2, s3, sx, tx, sy, ty, tz
#  Bad: t1, s2, t3, sz
# *** Good maps for P3 -> P2 extension:
[[HN(twist(P2), stab), HN(twist(P3), stab), HN(twist(s3(P2)), stab)] for twist in [s1,t1,s2,t2,s3,t3,sx,tx, sy, ty, sz, tz]]
# Good: s1, t1, t2, s3, sx, tx, sy, ty, sz, tz
#  Bad: s2, t3
# *** Good maps for P1 -> P2 -> P3 -> P3 extension:
[[HN(twist(z), stab), HN(twist(P3), stab), HN(twist(t3(z)), stab)] for twist in [s1,t1,s2,t2,s3,t3,sx,tx, sy, ty, sz, tz]]
# Bad: s3, sy, ty, tz
# *** Good maps for P1 -> P1 -> P2 -> P3 extension:
[[HN(twist(P1), stab), HN(twist(z), stab), HN(twist(s1(z)), stab)] for twist in [s1,t1,s2,t2,s3,t3,sx,tx, sy, ty, sz, tz]]
# Bad: t1, sz
# *** Guess for atoms:
def doesLower(o, twist, stab):
    return all([phase(h) < phase(o) for h in HN(twist(o),stab)])

#level4 = [[s1,s2,s3,s1]]
#level3 = [[s1,s2,s3]; [s1, s2, s1] ]
#level2 = [[s1,s2], [s1,sy], [s1, s3], [s2, s3]; [s1, s2], [s1, s1], [s2, s1]; ...]
#level1 = [s1, s2, s3, sx, sy, sz]
# *** Good aisles:
# <= phi is a good aisle for beta if beta preserves it.
def isGoodAisle(twist, stable, stab):
    # Here phi = phi(stable)
    for o in stab:
        if phase(o, stab) <= phase(stable, stab) and not all([phase(piece, stab) <= phase(stable, stab) for piece in HN(twist(o), stab)]):
            return False
        elif not all([phase(piece, stab) <= phase(stable, stab) for piece in HN(twist(o.shift(-1)), stab)]):
            return False
    return True

# <= phi is a better aisle for beta if beta strictly lowers it.
def isBetterAisle(twist, stable, stab):
    # Here phi = phi(stable)
    for o in stab:
        if phase(o, stab) <= phase(stable, stab) and not all([phase(piece, stab) < phase(stable, stab) for piece in HN(twist(o), stab)]):
            return False
        elif not all([phase(piece, stab) < phase(stable, stab) for piece in HN(twist(o.shift(-1)), stab)]):
            return False
    return True

twists = [s3,sy,sz,s1,sx,s2]
names = ['p3', 'y', 'z', 'p1', 'x', 'p2']

doubleSignatures = {(i,j):[isBetterAisle(composeAll([twists[i], twists[j]]), o, stab) for o in stab] for i in range(0,6) for j in range(0,6)}
tripleSignatures = {(i,j,k):[isBetterAisle(composeAll([twists[i], twists[j], twists[k]]), o, stab) for o in stab] for i in range(0,6) for j in range(0,6) for k in range(0,6)}
