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

