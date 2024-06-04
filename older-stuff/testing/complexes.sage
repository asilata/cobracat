load("../complexes.sage")

k = QQ
F = ProjectiveModuleOverField(k, 1)

# Test for checking a map.
P = ProjectiveComplex(k, {0: [F, F], 1:[F]}, {0:matrix(k, [[1], [-1]]).dict()}, {F: 'k'})
Map = {0: (3*identity_matrix(k,2)).dict(), 1:(3*identity_matrix(k,1)).dict()}
Trap = {0: (3*identity_matrix(k,2)).dict(), 1:(2*identity_matrix(k,1)).dict()}

# Test for direct sum.
Q = ProjectiveComplex(k, {-2: [F, F, F], -1:[F], 1: [F]}, {-2:matrix(k, [[2],[3],[1]]).dict()}, {F:'F'})
Q2 = Q.shift(2)
Qm2 = Q.shift(-2)

# Test for cone
C = cone(P, P, Map)

