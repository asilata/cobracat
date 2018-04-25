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
P = ProjectiveComplex(k, {0: [F, F], 1:[F]}, {0:matrix([[1, -1]]).dict()}, {F: 'k'})
Map = {0: (3*identity_matrix(2)).dict(), 1:(3*identity_matrix(1)).dict()}
Trap = {0: (3*identity_matrix(2)).dict(), 1:(2*identity_matrix(1)).dict()}

