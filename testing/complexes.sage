load("../complexes.sage")

k = QQ
R.<x> = PolynomialRing(k,1)
F = FreeModule(R, 1)

P = ProjectiveComplex(R)
P.addObject(0,F,'F')
P.addObject(1,F)
P.addObject(1,F)
P.addObject(-3,F)
P.addMap(0,0,0,x)
P.checkComplexity()
P.addMap(0, 0, 1, 1)
P.addObject(0, F)
P.addObject(0, F)
P.addObject(1, F)
P.addMap(0,0,1,2)
P.addMap(0,1,1,1)

# Test for checking a map.
F = FreeModule(k, 1)
P = ProjectiveComplex(k, {0: [F, F], 1:[F]}, {0:matrix(k, [[1], [-1]]).dict()}, {F: 'k'})
Map = {0: (3*identity_matrix(k,2)).dict(), 1:(3*identity_matrix(k,1)).dict()}
Trap = {0: (3*identity_matrix(k,2)).dict(), 1:(2*identity_matrix(k,1)).dict()}

# Test for direct sum.
Q = ProjectiveComplex(k, {-2: [F, F, F], -1:[F], 1: [F]}, {-2:matrix(k, [[2],[3],[1]]).dict()}, {F:'F'})
Q2 = Q.shift(2)
Qm2 = Q.shift(-2)

# Test for cone
C = cone(P, P, Map)

