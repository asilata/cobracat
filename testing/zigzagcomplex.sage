laffioad("../complexes.sage")
load("../zigzagalgebra.sage")
load("../zigzagmodules.sage")
load("../braidactions.sage")

a3 = make_test(a3graph)
R = a3['Z']
R.inject_variables()
F = ZigZagModule(R, 1, name = "P1")
G = ZigZagModule(R, 2, name = "P2")
H = ZigZagModule(R, 3, name = "P3")

P1 = ProjectiveComplex(R)
P1.addObject(0, F)

P2 = ProjectiveComplex(R)
P2.addObject(0, G)

P3 = ProjectiveComplex(R)
P3.addObject(0, H)

from functools import reduce
def composeAll(list_of_functions):
    return reduce (lambda x,y : lambda t : x(y(t)), list_of_functions, lambda x: x)


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

def u(C):
    return t2(s1(s2(C)))

def ui(C):
    return t2(t1(s2(C)))

def w(C):
    return s2(s3(t2(C)))

def wi(C):
    return s2(t3(t2(C)))

def v(C):
    return s1(s2(s3(t2(t1(C)))))

def vi(C):
    return s1(s2(t3(t2(t1(C)))))

def actOnPs(b):
    return b(P1),b(P2),b(P3)

def actOnQs(b):
    return b(P1),b(s1(P2)),b(s1(s2(P3)))

#Q1 = sigmaInverse(R, 2, P)
#Q2 = sigmaInverse(R, 2, Q1)
#Q2.minimize()
#Q3 = sigmaInverse(R, 2, Q2)

#R1 = sigma(R, 1, P)
#R2 = sigma(R, 1, R1)

#C = s1(s1(s1(s1(s2(t1(s2(t1(s2(P)))))))))


