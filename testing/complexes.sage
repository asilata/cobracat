load("../complexes.sage")

k = GF(101)
R.<x> = PolynomialRing(k,1)
F = FreeModule(R, 1)

P = ProjectiveComplex()
P.addObject(0,F,'F')
P.addObject(1,F)
P.addObject(1,F)
P.addObject(-3,F)
P.addMap(0,0,0,x)
