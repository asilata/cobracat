from sage.rings.noncommutative_ideals import Ideal_nc
from itertools import product
class ZigZagIdeal(Ideal_nc):
    def __init__(self, R):
        # Kill all length two paths with different sources and targets.
        # Equate all length two loops from each vertex.
        lentwopaths = filter(lambda x: x != 0, [R.prod(m) for m in product(R.arrows(), repeat = 2)])
        lenthreepaths = filter(lambda x: x!= 0, [R.prod(m) for m in product(R.arrows(), repeat = 3)])
        relations = lenthreepaths
        for e in R.idempotents():
            for f in R.idempotents():
                paths = filter(lambda x: x != 0, [R.prod([e,x,f]) for x in lentwopaths])
                if e == f:
                    relations.extend([x - y for (x,y) in zip(paths, paths[1:])])
                else:
                    relations.extend(paths)
        Ideal_nc.__init__(self, R, relations)
        
    def reduce(self,x):
        R = self.ring()
        # UNIMPLEMENTED
        # Remove any paths of length at least three, or any length two paths that are not loops.
        # Given a loop, find a canonical representative for it (probably best to fix a total order for the vertices beforehand).
        # return add([c*R(m) for m,c in x if len(m)<self._power],R(0))

# Standard test case        
d = {1:{2: 'a'}, 2:{1:'b', 3:'c'}, 3:{2:'d'}}
G = DiGraph(d)
A = G.path_semigroup().algebra(QQ)
