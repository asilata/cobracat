load("../complexes.sage")
load("../zigzagalgebra.sage")

a2 = make_test(a2graph)
R = a2['Z']
F = 'P1'
x = R.gens()[2]

P = ProjectiveComplex(R)
P.addObject(0,F)
P.addObject(1,F)
P.addObject(1,F)
P.addObject(-3,F)
P.addMap(0,0,0,3)
#P.checkComplexity()
#P.addMap(0, 0, 1, 1)
#P.addObject(0, F)
#P.addObject(0, F)
#P.addObject(1, F)
#P.addMap(0,0,1,2)
#P.addMap(0,1,1,1)
#P.checkComplexity()


