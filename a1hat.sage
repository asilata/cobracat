from itertools import product
from sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra_element import FiniteDimensionalAlgebraElement
from sage.modules.module import Module
from sage.structure.element import ModuleElement

class A1HatAlgebraElement(FiniteDimensionalAlgebraElement):
    def __init__(self, A, elt=None, check=True):
        FiniteDimensionalAlgebraElement.__init__(self, A = A, elt = elt, check = check)


class A1HatAlgebra(FiniteDimensionalAlgebra):
    '''
    The class for the A1-hat zigzag algebra equipped with a basis.
    '''
    Element = A1HatAlgebraElement
    
    def __init__(self, k): # k = base field, P = path semigroup of a quiver
        '''
        The ZigZag Algebra over base field k.
        P is the path semigroup of a quiver.
        '''
        a1hatgraph = DiGraph({1: {2: ['x1', 'y1']}, 2:{1:['x2', 'y2']}}, multiedges=True)
        P = a1hatgraph.path_semigroup()
        R = P.algebra(k)
        self._path_semigroup = P
        self._basis = list(R.idempotents()) + list(R.arrows()) + getLoops(R)
        table = [_getMatrix(R, self._basis, x) for x in self._basis]
        names = [str(x).replace('*','') for x in self._basis]
        super(A1HatAlgebra, self).__init__(k, table, names, category=Algebras(k).FiniteDimensional().Graded().WithBasis().Associative())

    def _repr_(self):
        return "A1 Hat zigzag algebra of {0} over {1}".format(self._path_semigroup.quiver(), self._base)

    @cached_method
    def isA1Hat(self):
        return True

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
        An element a such that a*b and b*a are both loops. Note that deg(a) + deg(b) = 2.
        More explicitly, if b is a loop at a vertex v, then a is the idempotent at v; if b is the idempotent at v then a is a loop at v; if b represents an an arrow, then a represents the reverse arrow. The element b must be in the basis.
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
        if len(sorted(set((degs)))) == 1:
            return degs[0]
        else:
            raise Exception("Element not homogeneous: " + str(a))

# Helper functions to get zigzag algebra elements and relations from a path algebra.
#@cached_method
def a1HatReduce(R,x):
    """
    Reduce an element x of R modulo the zigzag relations.
    """
    x1,x2,y1,y2 = R.arrows()[0], R.arrows()[1], R.arrows()[2], R.arrows()[3]
    monomials = x.sort_by_vertices()
    reduction = 0
    for m in monomials:
        z,v1,v2 =m[0],m[1],m[2]
        reduction = reduction + add([c*R(t) for t,c in z if len(t) < 2]) # Keep everything of length less than 2.
        if v1 == v2:
            for t,c in z:
                if len(t) != 2:
                    continue
                if t in [x1*y2, x2*y1, y1*x2, y2*x1]:
                    continue
                if t == y1*y2:
                    t = x1*x2
                if t == y2*y1:
                    t = x2*x1
                reduction = reduction + c*R(t)
    return reduction

@cached_method
def getLoops(R):
    """
    Return self-loops (of length two) in the path algebra R.
    """
    x1,x2,y1,y2 = R.arrows()[0], R.arrows()[1], R.arrows()[2], R.arrows()[3]
    return [R(x1*x2), R(x2*x1)]

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
    reduction = R(a1HatReduce(R,x))
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

class A1HatModuleElement(ModuleElement):
    def __init__(self, parent, x):
        self.x = x
        ModuleElement.__init__(self, parent = parent)

class A1HatModule(Module):
    '''
    Projective (left) modules over the A1-hat zigzag algebra
    '''

    Element = A1HatModuleElement
    
    def __init__(self, R, i, twist = 0, name="P", basis=None):
        '''
        The projective module P = R*e, where e is the i-th idempotent in R.idempotents()
        '''
        self._ring = R
        self._i = i
        self._idempotent = R.idempotents()[i-1]  # Vertices are conventionally 1,2,3,... but list elements are zero-indexed :(
        self._twist = twist
        self._name = name
        self._basis = basis
        Module.__init__(self, R, category=LeftModules(R))

    def copy(self):
        return A1HatModule(self._ring, self._i, self._twist, self._name)
        
    def twistBy(self, n = 1):
        '''
        Degree twist by n.
        '''
        return A1HatModule(self._ring, self._i, self._twist + n, self._name)

    def __repr__(self):
        return self._name + "<" + str(self._twist) + ">"

    def __str__(self):
        return self.__repr__()

    def idempotent(self):
        return self._idempotent

    def twist(self):
        return self._twist

    def name(self):
        return self._name

    @cached_method
    def hom(self, Q):
        e = self.idempotent()
        f = Q.idempotent()
        return filter(lambda x: x != 0, [e*x*f for x in self._ring.basis()] )

    # Here the module is supposed to be the left R-module Re, but the map is right multiplication by r.
    @cached_method
    def is_zero(self, r):
        '''
        Is right multiplication by r the zero map on self?
        '''
        return (self._idempotent * r == 0)

    @cached_method
    def basis(self):
        return self._basis
    
    @cached_method
    def is_invertible(self, r):
        '''
        Is right multiplication by r an invertible map on self?
        '''
        ir = self._idempotent * r
        if ir == 0:
            return False

        nonZeroCoeffs = [(x,y) for (x,y) in self._idempotent.monomial_coefficients().items() if y != 0]
        m,d = nonZeroCoeffs[0]
        c = ir.monomial_coefficients().get(m, 0)
        multiple = c/d

        return (ir == multiple*self._idempotent)

    @cached_method
    def invert(self, r):
        '''
        An element s in basering such that right multiplication by s is the inverse of right multiplication by r.
        Assumes that right multiplication by r is invertible.
        '''
        if not self.is_invertible(r):
            raise TypeError("Not invertible")
        else:
            ir = self._idempotent * r
            nonZeroCoeffs = [(x,y) for (x,y) in self._idempotent.monomial_coefficients().items() if y != 0]
            m,d = nonZeroCoeffs[0]
            c = ir.monomial_coefficients().get(m, 0)
            multiple = c/d
            return d/c * self._idempotent


def s(i, C):
    D = sigma(R, i, C)
    D.minimize()
    return D

def t(i, C):
    D = sigmaInverse(R, i, C)
    D.minimize()
    return D
        
def s1(C):
    return s(1,C)

def s2(C):
    return s(2,C)

def t1(C):
    return t(1,C)

def t2(C):
    return t(2,C)
