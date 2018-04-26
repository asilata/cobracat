load("../complexes.sage")
load("../zigzagalgebra.sage")
load("../zigzagmodules.sage")
load("../braidactions.sage")

a3 = make_test(a3graph)
R = a3['Z']
F = ZigZagModule(R, 1, name = "P1")

P = ProjectiveComplex(R)
P.addObject(0, F)
Q1 = sigmaInverseComplex(R, 2, P)
Q2 = sigmaInverseComplex(R, 2, Q1)
Q2.minimize()
Q3 = sigmaInverseComplex(R, 2, Q2)



