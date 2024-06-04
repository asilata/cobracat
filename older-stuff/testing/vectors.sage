# Loading this file will provide the rootplot function.
# Rootplot (a,b,c) plots the Z of the root (a,b,c).
# The stability condition is such that
# Z : alpha1 = (-1,0)
#   : alpha2 = (0, 1+e2)
#   : alpha3 = (1-e1, e1)
# These can be changed by changing the assignments at the beginning of the file.

e1 = 0.1
e2 = 4

a1 = vector([-1,0])
a2 = vector([0,1+e2])
a3 = vector([1-e1,e1])

M = matrix([a1,a2,a3]).transpose()

def rootplot(roots, main_root = vector([0,0,0]), names = None, **args):
    vectors = [M * vector(s) for s in roots]
    plots = [v.plot(**args) for v in vectors]
    if names == None:
        names = roots
    labels = [text(str(names[i]), vectors[i]) for i in range(0, len(roots))]
    return sum(plots) + sum(labels) + (M * vector(main_root)).plot(color='red')


def area_polygon(vecs):
    if len(vecs) < 3:
        return 0
    sum = 0
    for i in range(0,len(vecs)-1):
        v1,v2 = vecs[i],vecs[i+1]
        sum = sum + v1[0]*v2[1] - v1[1]*v2[0]
    sum = sum + vecs[-1][0]*vecs[0][1] - vecs[-1][1]*vecs[0][0]
    return sum/2

def subq_area(vecs):
    lst = [0 for i in range(0,len(vecs))]
    for i in range(0,len(vecs)):
        lst[i] = sum(vecs[0:i+1])
    lst = [(0,0)] + lst
    return area_polygon(lst)
