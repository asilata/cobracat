load("../complexes.sage")
load("../zigzagalgebra.sage")
load("../zigzagmodules.sage")
load("../braidactions.sage")

a3affgraph = DiGraph({1: {2: 'a1', 3: 'c2'}, 2:{1: 'a2', 3:'b1'}, 3:{2:'b2', 1:'c1'}})
a3aff = make_test(a3affgraph)
R = a3aff['Z']
R.inject_variables()

F1 = ZigZagModule(R, 1, name="P1")
F2 = ZigZagModule(R, 2, name="P2")
F3 = ZigZagModule(R, 3, name="P3")

P1 = ProjectiveComplex(R)
P1.addObject(0, F1)

P2 = ProjectiveComplex(R)
P2.addObject(0, F2)

P3 = ProjectiveComplex(R)
P3.addObject(0, F3)

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

u = composeAll([s1,s2,t1])
x = composeAll([t3,s1,s2,t1,t3])
y = composeAll([s2,t3,s1,s2,t1,s3,t2])


objects = [P1, P2, P3, s1(P2), s1(P3), s2(P3)]

def actOnPs(b):
    return b(P1),b(P2),b(P3)

def actOnAll(b):
    return b(P1),b(P2),b(P3),b(s1(P2)), b(s2(P3)), b(s1(P3))

def actOnQs(b):
    return b(s1(s3(s2(P1)))), b(s1(s3(s2(t1(P3))))), b(s1(s3(s2(t1(t3(P2)))))), b(s1(s3(s2(t1(t3(t2(P3))))))), b(s1(s3(s2(t1(t3(t2(t3(P1))))))))


