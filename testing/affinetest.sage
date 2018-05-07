load("../complexes.sage")
load("../zigzagalgebra.sage")
load("../zigzagmodules.sage")
load("../braidactions.sage")

a4affgraph = DiGraph({1: {2: 'a1', 4: 'd2'}, 2:{1: 'a2', 3:'b1'}, 3:{2:'b2', 4:'c1'}, 4:{3: 'c2', 1:'d1'}})
a4aff = make_test(a4affgraph)
R = a4aff['Z']
R.inject_variables()

F1 = ZigZagModule(R, 1, name="P1")
F2 = ZigZagModule(R, 2, name="P2")
F3 = ZigZagModule(R, 3, name="P3")
F4 = ZigZagModule(R, 4, name="P4")

P1 = ProjectiveComplex(R)
P1.addObject(0, F1)

P2 = ProjectiveComplex(R)
P2.addObject(0, F2)

P3 = ProjectiveComplex(R)
P3.addObject(0, F3)

P4 = ProjectiveComplex(R)
P4.addObject(0, F4)

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

def b1(C):
    return s4(t3(s2(s1(s3(C)))))

def b2(C):
    return s2(t1(s4(s3(s1(C)))))

def actOnPs(b):
    return b(P1),b(P2),b(P3),b(P4)

              
