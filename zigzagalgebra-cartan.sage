from itertools import product
from sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra_element import FiniteDimensionalAlgebraElement

class ZigZagAlgebraElement(FiniteDimensionalAlgebraElement):
    def __init__(self, A, elt=None, check=True):
        FiniteDimensionalAlgebraElement.__init__(self, A = A, elt = elt, check = check)
        

class ZigZagAlgebra(FiniteDimensionalAlgebra):
    r"""
    The zig-zag algebra over a field associated to a Cartan type.
    """
    Element = ZigZagAlgebraElement
    
    def __init__(self, ct, k):
        r'''
        Initialize ``self``.

        k is the base field, and ct is a CartanType.
        '''
        #R = P.algebra(k)
        #self._path_semigroup = P
        ct = CartanType(ct)
        self._cartan_type = ct
        if not self._cartan_type.is_simply_laced():
            raise ValueError("Given Cartan Type must be simply laced.")

        # This is currently a list of string names, for the idempotents, degree one edges, and loops of the zigzag algebra respectively.
        self._basis = _zz_basis(ct)

        # A list of matrices, where the ith matrix is the matrix of right multiplication
        # by the ith basis element
        table = [_zz_right_multiply(self._basis, b, k) for b in self._basis]
        names = [_zz_get_name(b) for b in self._basis]
        
        super(ZigZagAlgebra, self).__init__(k, table, names, category=Algebras(k).FiniteDimensional().Graded().WithBasis().Associative())

    def _repr_(self):
        print("Unimplemented.")
        return None

    def isA1Hat(self):
        print("Unimplemented.")
        return None

    # def pathsFromTo(self,e,f):
    #     '''
    #     Returns the number of paths in the algebra that go from idempotent e to idempotent f.
    #     Returns an empty list if e or f is not a primitive idempotent in the algebra.
    #     '''
    #     ids = self.idempotents()
    #     if e not in ids or f not in ids:
    #         return []
    #     paths = [e * x * f for x in self.basis()]
    #     return [p for p in paths if p != 0]
    
    # def idempotents(self):
    #     '''
    #     The list of idempotents of self.
    #     '''
    #     return self.basis()[0:len(self._path_semigroup.idempotents())]

    # def coeff(self, r, monomial):
    #     '''
    #     The coefficient of monomial in the expansion of r in the basis. 
    #     Monomial must be a basis element.
    #     '''
    #     i = list(self.basis()).index(monomial)
    #     return r.monomial_coefficients().get(i, 0)

    # def arrows(self):
    #     '''
    #     The list of elements of this algebra representing the arrows in the underlying (doubled) quiver.
    #     '''
    #     return self.basis()[len(self._path_semigroup.idempotents()):len(self._path_semigroup.idempotents())+len(self._path_semigroup.arrows())]

    # def loops(self):
    #     '''
    #     The list of elements of this algebra representing the loops. Note that there is a unique loop at every vertex.
    #     '''
    #     return self.basis()[len(self._path_semigroup.idempotents()) + len(self._path_semigroup.arrows()):]

    # def source(self, b):
    #     '''
    #     The idempotent corresponding to the source vertex of the arrow represented by b.
    #     '''
    #     if b not in self.basis():
    #         raise Exception("{0} is not a basis element.".format(b))
    #     for e in self.idempotents():
    #         if e * b != 0:
    #             return e
    #     raise Exception("Something went wrong: {0} does not seem to have a head.".format(b))

    # def target(self, b):
    #     '''
    #     The idempotent corresponding to the target vertex of the arrow represented by b.
    #     '''
    #     if b not in self.basis():
    #         raise Exception("{0} is not a basis element.".format(b))
    #     for e in self.idempotents():
    #         if b * e != 0:
    #             return e
    #     raise Exception("Something went wrong: {0} does not seem to have a head.".format(b))

    # def dualize(self, b):
    #     '''
    #     The element a such that a*b and b*a are both loops. Note that deg(a) + deg(b) = 2.
    #     More explicitly, if b is a loop at a vertex v, then a is the idempotent at v; if b is the idempotent at v then a is the loop at v; if b represents an an arrow, then a represents the reverse arrow. The element b must be in the basis.
    #     '''
    #     if b not in self.basis():
    #         raise Exception("{0} is not a basis element.".format(b))
    #     s,t = self.source(b),self.target(b)
    #     for c in self.basis():
    #         if b * t * c * s in self.loops():
    #             return c

    # def dualPairs(self, e, f):
    #     auxPairs = [(path, self.dualize(path)) for path in [x * e for x in self.basis()] if path != 0]
    #     return [(t,v*f) for (t,v) in auxPairs]
    
    # # Returns the degree of a homogeneous element.
    # @cached_method    
    # def deg(self,a):
    #     '''
    #     Degree of a, assuming a is homogeneous.
    #     '''
    #     mons = a.monomials()
    #     if mons == []:
    #         return -1
    #     degs = []
    #     for x in mons:
    #         if x in self.idempotents():
    #             degs = degs + [0]
    #         elif x in self.arrows():
    #             degs = degs + [1]
    #         else:
    #             degs = degs + [2]
    #     if len(set(degs)) == 1:
    #         return degs[0]
    #     else:
    #         raise Exception("Element not homogeneous: " + str(a))

# Helper functions to get zigzag algebra elements and relations from a path algebra.
# def _getTwoSteps(R):
#     """
#     Return a dictionary with keys all pairs of idempotents of R. The value of each pair of idempotents is the list all length 2 paths beginning and ending at the corresponding vertices.
#     """
#     twosteps = {}
#     len2paths = list(filter(lambda x: x != 0, [R.prod(m) for m in product(R.arrows(), repeat = 2)]))
#     for e in R.idempotents():
#         for f in R.idempotents():
#             twosteps[(e,f)] = list(filter(lambda x: x != 0, [R.prod([e,x,f]) for x in len2paths]))
#     return twosteps

# def zigZagReduce(R,x):
#     """
#     Reduce an element x of R modulo the zigzag relations.
#     """
#     monomials = x.sort_by_vertices()
#     reduction = 0
#     twosteps = _getTwoSteps(R)
#     for m in monomials:
#         z,v1,v2 =m[0],m[1],m[2]
#         reduction = reduction + add([c*R(m) for m,c in z if len(m) < 2]) # Keep everything of length less than 2.
#         if v1 == v2: # If we're looking at loops,
#             e = _getIdempotent(R, v1) # get the corresponding idempotent, and
#             l = twosteps[(e,e)][0] # select the first 2-loop in the previously chosen order.
#             reduction = reduction + add([c*l for m,c in z if len(m) == 2]) # Then add up copies of this chosen loop.
#     return reduction

# def getLoops(R):
#     """
#     Return self-loops (of length two) in the path algebra R.
#     """
#     loops = []
#     twosteps = _getTwoSteps(R)
#     for e in R.idempotents():
#         eloops = twosteps[(e,e)]
#         if eloops != []:
#                 loops.append(zigZagReduce(R,eloops[0]))
#     return loops

# def _getIdempotent(R,v):
#     """
#     Return the idempotent of R corresponding to the vertex v.
#     """
#     vertices = R.quiver().vertices()
#     idempotents = R.idempotents()
#     alist = [(x,y) for (x,y) in zip(vertices,idempotents) if x == v]
#     if len(alist) != 1:
#         raise Exception("Vertex does not correspond to a unique idempotent!")
#     else:
#         return alist[0][1]

# def _getCoefficients(R, basis, x):
#     '''
#     Return the coefficients of x wrt the given basis, after a zigzag reduction.
#     '''
#     reduction = R(zigZagReduce(R,x))
#     coeffDict = {R(k):v for k,v in reduction.monomial_coefficients().items()}
#     return [coeffDict.get(b,0) for b in basis]

# def _getMatrix(R, basis, x):
#     '''
#     Return the matrix of right multiplication by x in the given basis, modulo the zig zag relations.
#     '''
#     xMatrix = []
#     for y in basis:
#         coeffs = _getCoefficients(R, basis, y*x)
#         xMatrix.append(coeffs)
#     return matrix(xMatrix)

def _zz_idempotents(ct):
    r"""
    Given a Cartan Type, generate idempotent names for the zigzag algebra for each vertex.
    """
    return [(i,0) for i in ct.index_set()]

def _zz_edges(ct):
    r"""
    Given a Cartan Type, generate names for length-one paths in the zigzag algebra.

    In the Dynkin diagram of the cartan type, assume that edges are ordered from the
    lower labelled vertex to the higher labelled vertex, for example, 3 --> 4.
    In this case, there will be two generators corresponding to this edge, namely 'a_3_4' and 'b_4_3'.

    We assume that the given Cartan Type is simply laced.
    """
    dynkin_edges = ct.dynkin_diagram().edges()
    forward_zigzag_edges = [((i,j),1) for (i,j,_) in dynkin_edges if i < j]
    backward_zigzag_edges = [((i,j),1) for (i,j,_) in dynkin_edges if j < i]
    return forward_zigzag_edges + backward_zigzag_edges

def _zz_loops(ct):
    r"""
    Given a Cartan Type, generate loop names for the zigzag algebra for each vertex.
    """
    return [(i,2) for i in ct.index_set()]    

def _zz_basis(ct):
    r"""
    Given a Cartan Type, return a list of names of the basis elements of the
    corresponding zigzag algebra.
    """
    return _zz_idempotents(ct) + _zz_edges(ct) + _zz_loops(ct)

def _zz_deg(b):
    r"""
    Given a basis element of the form (i,n) for n in {0,2} or ((i,j),1), return its degree.
    """
    return b[1]

def _zz_source(b):
    r"""
    Given a basis element (i.e. a path of length at most 2), return its source vertex.
    """
    if _zz_deg(b) == 0 or _zz_deg(b) == 2:
        return b[0]
    else:
        return b[0][0]

def _zz_target(b):
    r"""
    Given a basis element (i.e. a path of length at most 2), return its target vertex.
    """
    if _zz_deg(b) == 0 or _zz_deg(b) == 2:
        return b[0]
    else:
        return b[0][1]

def _zz_get_name(b):
    r"""
    Given a basis element (i.e, a path of length at most 2), return its name.
    """
    if _zz_deg(b) == 0:
        return 'e' + str(_zz_source(b))
    elif _zz_deg(b) == 1:
        return 'a' + str(_zz_source(b)) + str(_zz_target(b))
    elif _zz_deg(b) == 2:
        return 'l' + str(_zz_source(b))
    else:
        raise ValueError("{} is not a valid basis element!".format(b))
                 
    
def _zz_right_multiply(basis, x, k):
    r"""
    Generate the multiplication table for right multiplication by a basis element x of the zigzag algebra.
    The element x must be in the basis.
    """

    # Initialize an empty transposed matrix, which we fill in.
    mult_matrix = Matrix(k, len(basis), len(basis))

    # Case 1: x is an idempotent
    if _zz_deg(x) == 0:
        for i in range(len(basis)):
            b = basis[i]
            if _zz_target(b) == _zz_source(x):
                j = basis.index(b)
                mult_matrix[i,j] = 1

    # Case 2: x is an edge (degree-one element) 
    elif _zz_deg(x) == 1:
        for i in range(len(basis)):
            b = basis[i]
            if _zz_deg(b) == 0 and _zz_target(b) == _zz_source(x):
                # Multiplication with idempotent at source of x is x.
                j = basis.index(b)
                mult_matrix[i,j] = 1
            elif (_zz_deg(b) == 1 and
                  _zz_target(b) == _zz_source(x) and 
                  _zz_source(b) == _zz_target(x)):
                # Multiplication with the opposite edge produces the loop at the target of x.
                # CHECK SIGN CONVENTION HERE.
                j = basis.index((_zz_target(x),2))
                if _zz_source(b) < _zz_target(b):
                    # b is a "forwards" arrow
                    mult_matrix[i,j] = 1
                else:
                    # b is a "backwards" arrow
                    mult_matrix[i,j] = -1                    
    # Case 3: x is a loop
    elif _zz_deg(x) == 2:
        for i in range(len(basis)):
            b = basis[i]
            if _zz_deg(b) > 0 or (_zz_target(b) != _zz_source(x)):
                # Multiplication with positive degree elements as well as
                # idempotents that don't match the endpoint of x is zero.
                pass
            else:
                # Multiplication with the idempotent with the same endpoint yields x.
                j = basis.index(x)
                mult_matrix[i,j] = 1
    else:
        raise ValueError("{} is not a basis element of the zigzag algebra!".format(x))
    return mult_matrix

