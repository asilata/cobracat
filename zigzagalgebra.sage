from sage.rings.noncommutative_ideals import Ideal_nc
from itertools import product
class ZigZagIdeal(Ideal_nc):
    def __init__(self, R):
        len2paths = filter(lambda x: x != 0, [R.prod(m) for m in product(R.arrows(), repeat = 2)])
        len3paths = filter(lambda x: x!= 0, [R.prod(m) for m in product(R.arrows(), repeat = 3)])
        self._twosteps = {} # For each pair of idempotents, record all length 2 paths beginning and ending at those vertices.
        for e in R.idempotents():
            for f in R.idempotents():
                self._twosteps[(e,f)] = filter(lambda x: x != 0, [R.prod([e,x,f]) for x in len2paths])

        relations = len3paths # Kill all paths of length three.
        for e in R.idempotents():
            for f in R.idempotents():
                eTof = self._twosteps[(e,f)]
                if e == f:
                    relations.extend([x - y for (x,y) in zip(eTof, eTof[1:])]) # Equate all 2-loops starting at a vertex.
                else:
                    relations.extend(eTof) # Kill all length 2 paths starting and ending at different vertices.
        Ideal_nc.__init__(self, R, relations)
        
    def reduce(self,x):
        R = self.ring()
        monomials = x.sort_by_vertices()
        reduction = 0
        for m in monomials:
            z,v1,v2 =m[0],m[1],m[2]
            reduction = reduction + add([c*R(m) for m,c in z if len(m) < 2]) # Keep everything of length less than 2.
            if v1 == v2: # If we're looking at loops,
                e = self._getIdempotent(v1) # get the corresponding idempotent, and
                l = self._twosteps[(e,e)][0] # select the first 2-loop in the previously chosen order.
                reduction = reduction + add([c*l for m,c in z if len(m) == 2]) # Then add up copies of this chosen loop.
        return reduction

    def loops(self):
        R = self.ring()
        loops = []
        for e in R.idempotents():
            eloops = self._twosteps[(e,e)]
            if eloops != []:
                loops.append(self.reduce(eloops[0]))
        return loops
                
    def _getIdempotent(self,v):
        vertices = self.ring().quiver().vertices()
        idempotents = self.ring().idempotents()
        alist = [(x,y) for (x,y) in zip(vertices,idempotents) if x == v]
        if len(alist) != 1:
            raise Exception("Vertex does not correspond to a unique idempotent!")
        else:
            return alist[0][1]
                

# Standard test case        
ds = {1:{2: 'a'}, 2:{1:'b', 3:'c'}, 3:{2:'d'}}
G = DiGraph(ds)
A = G.path_semigroup().algebra(QQ)

