# Loading this file will do the following main things:

# 1) Define R to be  the Zigzag algebra for the a3 quiver 

# 2) Define P1, P2, P3 to be the projective modules corresponding to the roots alpha1, alpha2, alpha3, considered as objects in the homotopy category of complexes (concentrated in degree 0).

# 3) Define s1, s2, s3 to be the positive twists on the homotopy category of complexes in P1, P2, P3.

# 4) Define t1, t2, t3 to be the negative twists on the homotopy category of complexes in P1, P2, P3.

# 5) Define heart as the list [P1, P2, P3]

# This should allow all computations without having to go into the details of the implementation.

# Common calculations:

# 1) s2(P1) returns [-1]: P2<-1> -> P1 :[0]
# The two numbers on either end represent homological the degree. The numbers in < > represent internal degree shifts.

# 2) s2(s1(P2)) returns [-1]: P1<-1> :[-1]

# 3) f = composeAll([s1,s2,t1,t2]) defines s as the transformation obtained by composing s1, s2, t1, t2. Having defined s in this way, it can be applied to any object as before. For example, f(P1) returns [-1]: P1<-2> -> P2<-1> :[0].

# 4) map(f, heart) returns a list obtained by applying f to all the elements of heart.

# Less common calculations

# 1) If K is a complex, then K.shift(n) returns the complex K[n].

# 2) K.show() shows a picture of K. The heart levels are color coded.

# 3) K.objects(i) returns the list of objects in K in homological degree (i).

# 4) K.maps(i) returns the list of maps from the objects at degree i to the objects at degree (i+1). The format is as follows: (a,b,m) indicates that there is a map from the a-th object in K.objects(i) (0 indexed) to the b-th object in K.objects(i+1) and the map is given by multiplication by the ring element m.


load("../complexes.sage")
load("../zigzagalgebra.sage")
load("../zigzagmodules.sage")
load("../braidactions.sage")
load("../HN.sage")

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

heart = [P1, P2, P3]

# def isUnimodal(p):
#     coeffs = p.coefficients()
#     if coeffs == []:
#         return True
#     increasing = True
#     current = coeffs[0]
#     peak = coeffs[0]
#     for c in coeffs[1:]:
#         if increasing:
#             if c[0] >= current[0]:
#                 current = c
#                 peak = c
#             else:
#                 increasing = False
#                 current = c
#         else:
#             if c[0] <= current[0]:
#                 current = c
#             else:
#                 print "Sequence not unimodal at power", current[1], "with coefficient", current[0]
#                 return False
#     print "Sequence unimodal with peak",peak[0], "at power",peak[1]
#     return True

# N = 100
# obDict = {}
# def constructOb():
#     """Populate obDict with (a,v) where v is the a'th computed point in the orbit. The structure is a binary tree where the left child of a node is the value of s1 and the right child of a node is the value of t2."""
#     c = 1
#     while c < N:
#         if c == 1:
#             obDict[1] = P1
#             c = c + 1
#             continue
#         z = obDict[c // 2]
#         if c % 2 == 0:
#             obDict[c] = s1(z)
#         else:
#             obDict[c] = t2(z)
#         c = c+1

# constructOb()
# vs = obDict.values()
# i = -1
# for x in vs:
#     i = i + 1
#     h1,h2 = hom(P1,x),hom(P2,x)
#     p1,p2 = h1.qPolynomial(), h2.qPolynomial()
#     if not isUnimodal(p1):
#         print "number of P2 not unimodal:",i, x,p1
#         break
#     if not isUnimodal(p2):
#         print "number of P1 not unimodal:",i, x,p2
#         break
    

