e1 = 0.1
e2 = 0.2

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

root_vec = (4,3,3)
root_seq = [(1,0,0),(1,0,1),(2,1,1),(2,1,2),(3,2,2),(3,2,3)]
proj = [M*vector(s) for s in root_seq]

#S = sum([(M * vector(s)).plot() for s in root_seq]) + sum([]) + (M * vector([4,3,3])).plot(color='red')

#mr = (2,1,2)
#rs = [(1,0,0),(1,0,1),(0,0,1)]

mr = (2,1,1)
rs = [(1,0,0),(1,0,1)]
