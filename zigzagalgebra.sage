from sage.rings.noncommutative_ideals import Ideal_nc
from itertools import product
from sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra_element import FiniteDimensionalAlgebraElement

class ZigZagAlgebraElement(FiniteDimensionalAlgebraElement):
    # def _div_(self, y):
    #     """
    #     If ``y`` is invertible, then return the quotient of ``self`` by ``y``.
    #     """
    #     if y.is_invertible():
    #         return x*y.inverse()
    #     else:
    #         raise ZeroDivisionError("Element {0} is not invertible".format(y))
            
    def __init__(self, A):
        super(FiniteDimensionalAlgebraElement, self).__init__(self, A)
    

class ZigZagAlgebra(FiniteDimensionalAlgebra):
    # Element = ZigZagAlgebraElement
    
    def __init__(self, k, P): # k = base field, P = path semigroup of a quiver
        R = P.algebra(k)
        I = ZigZagIdeal(R)
        self._path_semigroup = P
        self._basis = list(R.idempotents()) + list(R.arrows()) + I.loops()
        table = [_getMatrix(I, self._basis, x) for x in self._basis]
        names = [str(x).replace('*','') for x in self._basis]
        super(ZigZagAlgebra, self).__init__(k, table, names, category=Algebras(k).FiniteDimensional().WithBasis().Associative())

    def idempotents(self):
        return self.basis()[0:len(self._path_semigroup.idempotents())]

    def coeff(self, r, monomial):
        i = self.basis().index(monomial)
        return r.monomial_coefficients().get(i, 0)

    def arrows(self):
        return self.basis()[len(self._path_semigroup.idempotents()):len(self._path_semigroup.idempotents())+len(self._path_semigroup.arrows())]

    def loops(self):
        return self.basis()[len(self._path_semigroup.idempotents()) + len(self._path_semigroup.arrows()):]

    def source(self, b):
        if b not in self.basis():
            raise Exception("{0} is not a basis element.".format(b))
        for e in self.idempotents():
            if e * b != 0:
                return e
        raise Exception("Something went wrong: {0} does not seem to have a head.".format(b))

    def target(self, b):
        if b not in self.basis():
            raise Exception("{0} is not a basis element.".format(b))
        for e in self.idempotents():
            if b * e != 0:
                return e
        raise Exception("Something went wrong: {0} does not seem to have a head.".format(b))
    
    def dualize(self, b):
        if b not in self.basis():
            raise Exception("{0} is not a basis element.".format(b))
        s,t = self.source(b),self.target(b)
        for c in self.basis():
            if b * t * c * s in self.loops():
                return c
    
    # Returns the degree of a homogeneous element.
    def deg(self,a):
        mons = a.monomials()
        if mons == []:
            return -1
        degs = []
        for x in mons:
            if x in self.idempotents():
                degs = degs + [0]
            elif x in self.arrows():
                degs = degs + [1]
            else:
                degs = degs + [2]
        if len(uniq(degs)) == 1:
            return degs[0]
        else:
            raise Exception("Element not homogeneous: " + str(a))

# Returns the coefficients of x wrt the given basis, after reducing modulo the ideal I.
def _getCoefficients(I, basis, x):
    R = I.ring()
    coeffDict = {R(k):v for k,v in R(I.reduce(x)).monomial_coefficients().items()}
    return [coeffDict.get(b,0) for b in basis]

# Returns the matrix of right multiplication by x in the given basis, modulo the ideal I.
def _getMatrix(I, basis, x):
    R = I.ring()
    xMatrix = []
    for y in basis:
        xMatrix.append(_getCoefficients(I, basis, y*x))
    return matrix(xMatrix)
    
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

                

# Tests
# Test constructor
def make_test(graph, k=QQ):
    test = {}
    test['A'] = graph.path_semigroup().algebra(k)
    test['I'] = ZigZagIdeal(test['A'])
    test['Z'] = ZigZagAlgebra(k,graph.path_semigroup())
    return test

# Some standard graphs
a2graph = DiGraph({1: {2: 'a'}, 2:{1:'b'}})
a3graph = DiGraph({1:{2: 'a'}, 2:{1:'b', 3:'c'}, 3:{2:'d'}})
d4graph = DiGraph({1:{2:'a', 3:'b', 4:'c'}, 2:{1:'d'}, 3:{1:'e'}, 4:{1:'f'}})



