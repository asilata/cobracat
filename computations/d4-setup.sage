load("../complexes.sage")
load("../zigzagalgebra.sage")
load("../zigzagmodules.sage")
load("../braidactions.sage")

d4 = make_test(d4graph)
R = d4['Z']
R.inject_variables()
F = ZigZagModule(R, 1, name = "P1")
G = ZigZagModule(R, 2, name = "P2")
H = ZigZagModule(R, 3, name = "P3")
J = ZigZagModule(R, 4, name = "P4")

P1 = ProjectiveComplex(R)
P1.addObject(0, F)

P2 = ProjectiveComplex(R)
P2.addObject(0, G)

P3 = ProjectiveComplex(R)
P3.addObject(0, H)

P4 = ProjectiveComplex(R)
P4.addObject(0, J)


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

def s4(C):
    return s(4, C)

def t1(C):
    return t(1, C)

def t2(C):
    return t(2, C)

def t3(C):
    return t(3, C)

def t4(C):
    return t(4, C)

R = RootSystem(CartanType(['D',4]).relabel({1:1,2:4,4:2,3:3}))
R.root_lattice().dynkin_diagram()

central_charge = matrix([[-3,1],[-2,1], [-1,1],[7,0]])

change_of_basis = matrix([vector(R.ambient_space()(v)) for v in R.root_lattice().simple_roots()])

central_change_euclidean = change_of_basis.inverse() * central_charge
