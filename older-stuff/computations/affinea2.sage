# Loading this file will do the following main things:

# 1) Define R to be  the Zigzag algebra for the affinea2 quiver (rank 3).

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

heart = [P1, P2, P3]

def generateNextLevel(multipliers, given):
    '''
    Given a list of braids as "given", generate the next list by multiplying everything in given by everything in multipliers
    on the left.
    The givens is a list of list of keys in "generators", and multipliers is a list of keys in "generators"
    '''
    output = []
    for b in multipliers:
        output = output + [[b] + g for g in given]
    return output

def braidsGetOrAdd(n, dct, gens):
    '''
    Get keyval n from dct. If there is no such key, add braids of length n by getting/creating braids of length (n-1) and multiplying each one by each gen in gens.
    '''
    if n < 1:
        return None
    if n in dct.keys():
        return dct[n]
    elif n == 1:
        new = [[g] for g in gens]
    else:
        prev = braidsGetOrAdd(n-1, dct, gens)
        new = generateNextLevel(gens, prev)

    dct[n] = new
    return new

def allBraidsGetOrAdd(n):
    return braidsGetOrAdd(n, allBraidsByLen, generators.keys())

def positiveBraidsGetOrAdd(n):
    return braidsGetOrAdd(n, positiveMonoid, positiveGens)

def negativeBraidsGetOrAdd(n):
    return braidsGetOrAdd(n, negativeMonoid, negativeGens)
    
# def is_mixed(b):
#     if len(b) <=1:
#         return False
#     if b[0] in positiveGens:
#         fst, rest = positiveGens, negativeGens
#     else:
#         fst, rest = negativeGens, positiveGens

#     i = 0
#     for j in range(0,len(b) + 1):
#         if b[j] in fst:
#             i = j
#         else:
#             break
#     if j == len(b):
#         return True
#     else:
#         for 


# Test the definition of is_nonexpanding
# DEBUG = False

# if DEBUG:
#     assert is_nonexpanding(t3, t3(s2(P1))) == False
#     assert is_nonexpanding(s3, t3(s2(P1))) == True
#     assert is_nonexpanding(t2, t3(s2(P1))) == True
#     assert is_nonexpanding(s2, t3(s2(P1))) == False


