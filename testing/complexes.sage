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



