from functools import cached_property
from itertools import product
from sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra_element import FiniteDimensionalAlgebraElement

def _zz_idempotents(ct):
    r"""
    Given a Cartan Type, generate internal representations for the idempotents of
    the zigzag algebra for each vertex.

    We represent each idempotent as (i,0) where i is the index of the vertex, and 0
    represents that it is an element of degree 0.
    """
    return [(i,0) for i in ct.index_set()]

def _zz_arrows(ct):
    r"""
    Given a Cartan Type, generate internal representations for length-one paths in
    the zigzag algebra.

    In the Dynkin diagram of the cartan type, assume that edges are ordered from the
    lower labelled vertex to the higher labelled vertex, for example, 3 -> 4.
    In this case, there will be two generators corresponding to this edge, one 'forwards'
    and the other 'backwards'.

    We represent each arrow as ((i,j),0) where i is the source, j the target, and 1
    represents that it is an element of degree 1.

    We assume (for now) that the given Cartan Type is simply laced.
    """
    dynkin_edges = ct.dynkin_diagram().edges()
    forward_zigzag_arrows = [((i,j),1) for (i,j,_) in dynkin_edges if i < j]
    backward_zigzag_arrows = [((i,j),1) for (i,j,_) in dynkin_edges if j < i]
    return forward_zigzag_arrows + backward_zigzag_arrows

def _zz_loops(ct):
    r"""
    Given a Cartan Type, generate internal representations for the loops in
    the zigzag algebra for each vertex.

    We represent each loop as (i,2) where i is the source/target, and 2
    represents that it is an element of degree 2.
    """
    return [(i,2) for i in ct.index_set()]    

def _zz_basis(ct):
    r"""
    Given a Cartan Type, return a list of internal representations of the
    basis elements of the corresponding zigzag algebra.
    """
    return _zz_idempotents(ct) + _zz_arrows(ct) + _zz_loops(ct)

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
    Given a basis element (i.e, a path of length at most 2), construct its name.
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

        # For the moment we must input either a genuine CartanType, or a string shorthand
        # such as "A5" or "E8".
        # For some reason Sage does not seem to accept list inputs such as ['A', 5,1].
        self._cartan_type = CartanType(ct)
        self._base_ring = k

        if not self._cartan_type.is_simply_laced():
            # For the moment we only implement simply laced Cartan types.
            raise ValueError("Given Cartan Type must be simply laced.")

        # An internal representation for the basis of `self`, as specified in the helper methods `_zz_basis()`defined previously.
        self._internal_basis = _zz_basis(self._cartan_type)

        # A list of matrices, where the ith matrix is the matrix of right multiplication
        # by the ith basis element
        table = [_zz_right_multiply(self._internal_basis, b, k) for b in self._internal_basis]
        names = [_zz_get_name(b) for b in self._internal_basis]
        
        super(ZigZagAlgebra, self).__init__(k, table, names, category=Algebras(k).FiniteDimensional().Graded().WithBasis().Associative())

    def _repr_(self):
        return "Zig-zag algebra of {0} over {1}".format(self._cartan_type, self._base_ring)

    @property
    def _basis_correspondence(self):
        """
        Returns a zipped list of the internal representations of the basis elements of `self`
        with the algebra representations of the basis elements of `self`.
        """
        return zip(self._internal_basis, self.basis())

    def _to_internal_repr(self,b):
        """
        Returns the internal representation of a basis element `b` of `self`.

        Assume that `b` is in `self.basis()`.
        """
        if b not in self.basis():
            raise ValueError("{} is not a basis element of {}!".format(b,self))

        # Assume that _basis_correspondence is a 1-1 correspondence, so return the first element
        # of the first tuple whose second element is b.
        return [x for (x,y) in self._basis_correspondence if y == b][0]

    def _from_internal_repr(self,bi):
        """
        Returns the basis element corresponding to the internally represented basis element `bi`.

        Assume that `bi` is in self._internal_basis
        """
        if bi not in self._internal_basis:
            raise ValueError("{} is not the internal representation of any basis element of {}!".format(bi,self))

        # Assume that _basis_correspondence is a 1-1 correspondence, so return the second element
        # of the first tuple whose first element is bi.
        return [y for (x,y) in self._basis_correspondence if x == bi][0]

    def deg(self, p):
        """
        Returns the degree of a homogeneous element `p` of `self`.

        Assume that `p` is homogeneous.
        """
        mons = p.monomials()
        if mons == []:
            return -1
        degs = [_zz_deg(self._to_internal_repr(m)) for m in mons]
        for d in degs:
            if d != degs[0]:
                raise ValueError("{} is not a homogeneous element!".format(p))
        return degs[0]

    def source(self,b):
        """
        Returns the source vertex of the basis element `b` of `self`.

        Assume that `b` is in `self.basis()`.
        """
        if b not in self.basis():
            raise ValueError("{} is not a basis element of {}!".format(b,self))

        return _zz_source(self._to_internal_repr(b))

    def target(self,b):
        """
        Returns the target vertex of the basis element `b` of `self`.

        Assume that `b` is in `self.basis()`.
        """
        if b not in self.basis():
            raise ValueError("{} is not a basis element of {}!".format(b,self))
        
        return _zz_target(self._to_internal_repr(b))

    @cached_property
    def vertices(self):
        """
        Returns the index set of the Cartan Type of `self`. That is, the list of vertices of the Dynkin diagram of `self`.
        """
        return self._cartan_type.index_set()

    @cached_property
    def idempotents(self):
        '''
        The list of idempotents of `self`.
        '''
        return [x for x in self.basis() if self.deg(x) == 0]

    @cached_property
    def arrows(self):
        '''
        The list of arrows, or length-1 paths, in `self`.
        '''
        return [x for x in self.basis() if self.deg(x) == 1]        

    @cached_property
    def loops(self):
        '''
        The list of elements of this algebra representing the loops. Note that there is a unique loop at every vertex.
        '''
        return [x for x in self.basis() if self.deg(x) == 2]

    def idempotent_by_vertex(self, v):
        """
        Returns the idempotent of `self` corresponding to vertex v.
        """
        if v not in self.vertices:
            raise ValueError("{} is not a vertex of {}!".format(v,self))
        # Return the first (and only) idempotent that has v as a source vertex.
        return [e for e in self.idempotents if self.source(e) == v][0]

    # def isA1Hat(self):
    #     print("Unimplemented.")
    #     return None

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
    

    # def coeff(self, r, monomial):
    #     '''
    #     The coefficient of monomial in the expansion of r in the basis. 
    #     Monomial must be a basis element.
    #     '''
    #     i = list(self.basis()).index(monomial)
    #     return r.monomial_coefficients().get(i, 0)



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
    

