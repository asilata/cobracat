from itertools import product
from sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra_element import FiniteDimensionalAlgebraElement

class ZigZagAlgebraElement(FiniteDimensionalAlgebraElement):
    def __init__(self, A, elt=None, check=True):
        FiniteDimensionalAlgebraElement.__init__(self, A = A, elt = elt, check = check)
        

class ZigZagAlgebra(FiniteDimensionalAlgebra):
    '''
    The class for the ZigZag Algebra equipped with a basis.
    '''
    Element = ZigZagAlgebraElement
    
    def __init__(self, k, P): # k = base field, P = path semigroup of a quiver
        '''
        The ZigZag Algebra over base field k.
        P is the path semigroup of a quiver.
        '''
        R = P.algebra(k)
        self._path_semigroup = P
        self._basis = list(R.idempotents()) + list(R.arrows()) + getLoops(R)
        table = [_getMatrix(R, self._basis, x) for x in self._basis]
        names = [str(x).replace('*','') for x in self._basis]
        super(ZigZagAlgebra, self).__init__(k, table, names, category=Algebras(k).FiniteDimensional().Graded().WithBasis().Associative())

    def _repr_(self):
        return "Zig-zag algebra of {0} over {1}".format(self._path_semigroup.quiver(), self._base)

    @cached_method
    def isA1Hat(self):
        return False

    @cached_method
    def pathsFromTo(self,e,f):
        '''
        Returns the number of paths in the algebra that go from idempotent e to idempotent f.
        Returns an empty list if e or f is not a primitive idempotent in the algebra.
        '''
        ids = self.idempotents()
        if e not in ids or f not in ids:
            return []
        paths = [e * x * f for x in self.basis()]
        return [p for p in paths if p != 0]
    
    @cached_method
    def idempotents(self):
        '''
        The list of idempotents of self.
        '''
        return self.basis()[0:len(self._path_semigroup.idempotents())]

    @cached_method
    def coeff(self, r, monomial):
        '''
        The coefficient of monomial in the expansion of r in the basis. 
        Monomial must be a basis element.
        '''
        i = list(self.basis()).index(monomial)
        return r.monomial_coefficients().get(i, 0)

    @cached_method
    def arrows(self):
        '''
        The list of elements of this algebra representing the arrows in the underlying (doubled) quiver.
        '''
        return self.basis()[len(self._path_semigroup.idempotents()):len(self._path_semigroup.idempotents())+len(self._path_semigroup.arrows())]

    @cached_method
    def loops(self):
        '''
        The list of elements of this algebra representing the loops. Note that there is a unique loop at every vertex.
        '''
        return self.basis()[len(self._path_semigroup.idempotents()) + len(self._path_semigroup.arrows()):]

    @cached_method
    def source(self, b):
        '''
        The idempotent corresponding to the source vertex of the arrow represented by b.
        '''
        if b not in self.basis():
            raise Exception("{0} is not a basis element.".format(b))
        for e in self.idempotents():
            if e * b != 0:
                return e
        raise Exception("Something went wrong: {0} does not seem to have a head.".format(b))

    @cached_method
    def target(self, b):
        '''
        The idempotent corresponding to the target vertex of the arrow represented by b.
        '''
        if b not in self.basis():
            raise Exception("{0} is not a basis element.".format(b))
        for e in self.idempotents():
            if b * e != 0:
                return e
        raise Exception("Something went wrong: {0} does not seem to have a head.".format(b))

    @cached_method
    def dualize(self, b):
        '''
        The element a such that a*b and b*a are both loops. Note that deg(a) + deg(b) = 2.
        More explicitly, if b is a loop at a vertex v, then a is the idempotent at v; if b is the idempotent at v then a is the loop at v; if b represents an an arrow, then a represents the reverse arrow. The element b must be in the basis.
        '''
        if b not in self.basis():
            raise Exception("{0} is not a basis element.".format(b))
        s,t = self.source(b),self.target(b)
        for c in self.basis():
            if b * t * c * s in self.loops():
                return c

    @cached_method
    def dualPairs(self, e, f):
        auxPairs = [(path, self.dualize(path)) for path in [x * e for x in self.basis()] if path != 0]
        return [(t,v*f) for (t,v) in auxPairs]
    
    # Returns the degree of a homogeneous element.
    @cached_method    
    def deg(self,a):
        '''
        Degree of a, assuming a is homogeneous.
        '''
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

# Helper functions to get zigzag algebra elements and relations from a path algebra.
@cached_method
def _getTwoSteps(R):
    """
    Return a dictionary with keys all pairs of idempotents of R. The value of each pair of idempotents is the list all length 2 paths beginning and ending at the corresponding vertices.
    """
    twosteps = {}
    len2paths = list(filter(lambda x: x != 0, [R.prod(m) for m in product(R.arrows(), repeat = 2)]))
    for e in R.idempotents():
        for f in R.idempotents():
            twosteps[(e,f)] = list(filter(lambda x: x != 0, [R.prod([e,x,f]) for x in len2paths]))
    return twosteps

@cached_method
def zigZagReduce(R,x):
    """
    Reduce an element x of R modulo the zigzag relations.
    """
    monomials = x.sort_by_vertices()
    reduction = 0
    twosteps = _getTwoSteps(R)
    for m in monomials:
        z,v1,v2 =m[0],m[1],m[2]
        reduction = reduction + add([c*R(m) for m,c in z if len(m) < 2]) # Keep everything of length less than 2.
        if v1 == v2: # If we're looking at loops,
            e = _getIdempotent(R, v1) # get the corresponding idempotent, and
            l = twosteps[(e,e)][0] # select the first 2-loop in the previously chosen order.
            reduction = reduction + add([c*l for m,c in z if len(m) == 2]) # Then add up copies of this chosen loop.
    return reduction

@cached_method
def getLoops(R):
    """
    Return self-loops (of length two) in the path algebra R.
    """
    loops = []
    twosteps = _getTwoSteps(R)
    for e in R.idempotents():
        eloops = twosteps[(e,e)]
        if eloops != []:
                loops.append(zigZagReduce(R,eloops[0]))
    return loops

@cached_method
def _getIdempotent(R,v):
    """
    Return the idempotent of R corresponding to the vertex v.
    """
    vertices = R.quiver().vertices()
    idempotents = R.idempotents()
    alist = [(x,y) for (x,y) in zip(vertices,idempotents) if x == v]
    if len(alist) != 1:
        raise Exception("Vertex does not correspond to a unique idempotent!")
    else:
        return alist[0][1]

def _getCoefficients(R, basis, x):
    '''
    Return the coefficients of x wrt the given basis, after a zigzag reduction.
    '''
    reduction = R(zigZagReduce(R,x))
    coeffDict = {R(k):v for k,v in reduction.monomial_coefficients().items()}
    return [coeffDict.get(b,0) for b in basis]

def _getMatrix(R, basis, x):
    '''
    Return the matrix of right multiplication by x in the given basis, modulo the zig zag relations.
    '''
    xMatrix = []
    for y in basis:
        coeffs = _getCoefficients(R, basis, y*x)
        xMatrix.append(coeffs)
    return matrix(xMatrix)

# Tests
# Test constructor
def make_test(graph, k=QQ):
    test = {}
    test['A'] = graph.path_semigroup().algebra(k)
    test['Z'] = ZigZagAlgebra(k,graph.path_semigroup())
    return test

# Some standard graphs
a2graph = DiGraph({1: {2: 'a'}, 2:{1:'b'}})
a3graph = DiGraph({1:{2: 'a'}, 2:{1:'b', 3:'c'}, 3:{2:'d'}})
d4graph = DiGraph({4:{1:'a', 2:'b', 3:'c'}, 2:{4:'d'}, 3:{4:'e'}, 1:{4:'f'}})



