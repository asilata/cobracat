# Loading this file will provide the rootplot function.
# Rootplot [list of (a,b,c)] plots the central charge of the roots (a,b,c).
# The stability condition is such that
# Z : alpha1 = (-1,0)
#   : alpha2 = (0, 1+e2)
#   : alpha3 = (1-e1, e1)
# These can be changed by changing the assignments at the beginning of the file.

# Useful features of plots: 
# Plots can be assigned to variables.
# P = rootplot([(1,0,0),(0,1,1), (1,1,0)])
# will not show anything, but it will save the graph in P.
# P.show() 
# will then show the graph.
# Plots can be added.
# Q = P + rootplot([(0,0,1)])
# will create a new plot Q containing P and the additional vector corresponding to (0,0,1).
# Again, Q.show() will show it.

e1 = 0.5
e2 = 0

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

