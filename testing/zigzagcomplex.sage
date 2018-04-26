load("../complexes.sage")
load("../zigzagalgebra.sage")
load("../zigzagmodules.sage")
load("../braidactions.sage")

a3 = make_test(a3graph)
R = a3['Z']
F = ZigZagModule(R, 1, name = "P1")

P = ProjectiveComplex(R)
P.addObject(0, F)

def s1(C):
    D = sigma(R, 1, C)
    D.minimize()
    return D

def s2(C):
    D = sigma(R, 2, C)
    D.minimize()
    return D

def t1(C):
    D = sigmaInverse(R, 1, C)
    D.minimize()
    return D

def t2(C):
    D = sigmaInverse(R, 2, C)
    D.minimize()
    return D

Q1 = sigmaInverse(R, 2, P)
Q2 = sigmaInverse(R, 2, Q1)
Q2.minimize()
Q3 = sigmaInverse(R, 2, Q2)

R1 = sigma(R, 1, P)
R2 = sigma(R, 1, R1)




