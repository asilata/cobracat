from sage.rings.noncommutative_ideals import Ideal_nc
from itertools import product
class ZigZagIdeal(Ideal_nc):
    def __init__(self, R):
        # Kill all length two paths with different sources and targets.
        # Equate all length two loops from each vertex.
        lentwopaths = filter(lambda x: x != 0, [R.prod(m) for m in product(R.arrows(), repeat = 2)])
        self._twoloops = {}
        for e in R.idempotents():
            self._twoloops[e] = filter(lambda x: x != 0, [R.prod([e,x,e]) for x in lentwopaths])
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
        monomials = x.sort_by_vertices()
        reduction = 0
        for m in monomials:
            z,v1,v2 =m[0],m[1],m[2]
            reduction = reduction + add([c*R(m) for m,c in z if len(m) < 2])
            if v1 == v2:
                e = self.getIdempotent(v1)
                l = self._twoloops[e][0] # First loop in chosen order
                reduction = reduction + add([c*l for m,c in z if len(m) == 2])
        return reduction
                
    def getIdempotent(self,v):
        vertices = self.ring().quiver().vertices()
        idempotents = self.ring().idempotents()
        correspondence = [(x,y) for (x,y) in zip(vertices,idempotents) if x == v]
        if len(correspondence) != 1:
            raise Exception("Vertex does not correspond to an idempotent!")
        else:
            return correspondence[0][1]
                


# Standard test case        
d = {1:{2: 'a'}, 2:{1:'b', 3:'c'}, 3:{2:'d'}}
G = DiGraph(d)
A = G.path_semigroup().algebra(QQ)
