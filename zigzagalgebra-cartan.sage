r"""
Zig-zag algebras

Zig-zag algebras are (graded) finite-dimensional quotients of path algebras of quivers. Usually we consider zig-zag algebras of quivers without self-loops or multiple edges. In that case, the zig-zag algebra is a quotient of the path algebra of the doubled quiver by the following relations:
- all length-three paths are sent to zero,
- all length-two paths whose source is different from the target are sent to zero, and
- all length-two paths that start and end at the same vertex are set to be equal to each other up to sign. We call these length-two paths `loops`. The sign in the relations is determined by whether the first edge in the loop is an original edge in the quiver or a doubled edge.

There is a natural grading by path length as the relations are homogeneous. A basis for this algebra as a vector space is given by:
- length-zero paths at each vertex (the `idempotents`),
- length-one paths between two distinct vertices (the `arrows`), and
- for each vertex, any one back-and-forth path of length two starting and ending at that vertex (the `loops`).

AUTHORS:

- Asilata Bapat (2023-08-23): initial version

- Anand Deopurkar (2023-08-23): initial version

"""

# ****************************************************************************
#       Copyright (C) 2023 Asilata Bapat <asilata@alum.mit.edu> 
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from functools import cached_property
from itertools import product
from sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra_element import FiniteDimensionalAlgebraElement

def _zz_basis(ct):
    r"""
    Return a list of basis elements of the zigzag algebra of CartanType `ct`.
    
    The standard basis of the zigzag algebra consists of the primitive
    idempotents at each vertex, arrows between two vertices, and loop
    maps at each vertex.  Each path is represented as a tuple
    `(i,j,d)`, where `i` is its source vertex, `j` is its target
    vertex, and `d` is its degree.
    
    We assume (for now) that the given Cartan Type is simply laced.
    
    INPUT:

    - `ct` -- a CartanType

    OUTPUT:

    A list of paths that form the standard basis of the zigzag algebra. 
    """
    idempotents = [(i,i,0) for i in ct.index_set()]
    arrows = [(i,j,1) for (i,j,_) in ct.dynkin_diagram().edges()]
    loops = [(i,i,2) for i in ct.index_set()]
    
    return idempotents + arrows + loops

def _zz_source(b):
    r"""
    Return the source vertex of the given basis element.
    
    INPUT:

    - `b` -- the internal representation of a basis element of the zigzag algebra of CartanType `ct`.

    OUTPUT:

    The source vertex of the given basis element. Given (i,j,d), return i.
    """
    return b[0]

def _zz_target(b):
    r"""
    Return the target vertex of the given basis element.
    
    INPUT:

    - `b` -- the internal representation of a basis element of the zigzag algebra of CartanType `ct`.

    OUTPUT:

    The target vertex of the given basis element. Given (i,j,d), return j.
    """
    return b[1]

def _zz_deg(b):
    r"""
    Return the degree of a given basis element.

    INPUT:

    - `b` -- the internal representation of a  basis element of the zigzag algebra of CartanType `ct`.

    OUTPUT:

    The degree of the given basis element. Given (i,j,d), return d.
    """
    return b[2]

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

def _zz_right_multiply(basis, x, k=QQ):
    r"""
    Return the multiplication table for right multiplication by the basis element `x`.
    
    Generate the multiplication table for right multiplication of the elements in `basis` by a basis element `x` of the zigzag algebra.
    The element `x` must be in `basis`.
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
        Create a zig-zag algebra.

        INPUT:

        - `ct` -- CartanType
        - `k` -- Field

        OUTPUT:

        The zig-zag algebra of type `ct` over the field `k`.  The zig-zag algebra is
        `Associative`, `Graded`, `FiniteDimensional`, `Algebras(k)`, and `WithBasis`.
        
        '''

        # For the moment we must input either a genuine CartanType, or a string shorthand
        # such as "A5" or "E8".
        # For some reason Sage does not seem to accept list inputs such as ['A', 5,1].
        self.cartan_type = CartanType(ct)
        if not self.cartan_type.is_simply_laced():
            # For the moment we only implement simply laced Cartan types.
            raise ValueError("Given Cartan Type must be simply laced.")
        
        self.vertices = self.cartan_type.index_set()
        self._base_ring = k

        # An internal representation for the basis of `self`, as specified in the helper methods `_zz_basis()`defined previously.
        self._internal_basis = _zz_basis(self.cartan_type)

        # A list of matrices, where the ith matrix is the matrix of right multiplication
        # by the ith basis element
        table = [_zz_right_multiply(self._internal_basis, b, k) for b in self._internal_basis]
        names = [_zz_get_name(b) for b in self._internal_basis]
        
        super(ZigZagAlgebra, self).__init__(k, table, names, category=Algebras(k).FiniteDimensional().Graded().WithBasis().Associative())

    def _repr_(self):
        return "Zig-zag algebra of {0} over {1}".format(self.cartan_type, self._base_ring)

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
        Return the degree of the given homogeneous element.

        INPUT:

        - `p` -- A homogeneous `k`-linear combination of basis elements of `self`.

        OUTPUT:

        The degree of `p`.  Raise `ValueError` if `p` is not homogeneous.

        EXAMPLES:

           sage: A = ZigZagAlgebra("A2", QQ)
           sage: A.deg(A.basis()[0])
           0
           sage: A.deg(A.basis()[-1])
           2
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
        Returns the idempotent of `self` corresponding to vertex `v`.
        """
        if v not in self.vertices:
            raise ValueError("{} is not a vertex of {}!".format(v,self))
        # Return the first (and only) idempotent that has v as a source vertex.
        return [e for e in self.idempotents if self.source(e) == v][0]

    def calabi_yau_dual(self, b):
        r"""
        Return the unique Calabi--Yau dual basis element of `b`.

        INPUT:

        - `b`  -- a basis element of `self`.

        OUTPUT:

        - The unique basis element `a` of `self` such that `a * b` and `b * a` are both plus or minus of the loops at the corresponding vertices.
        
          `ValueError` if `b` is not in `self.basis()`.
        """
        if b not in self.basis():
            raise ValueError("{} is not a basis element of {}.".format(b,Z))
        s,t = self.source(b),self.target(b)
        for c in self.basis():
            bc = b * t * c * s
            if bc in self.loops() or -bc in self.loops():
                return c

    @cached_method
    def paths_from_to(self,e,f):
        r"""
        Returns a list of basis elements in the algebra that are paths from the source/target of `e` to the source/target of `f`.
        We assume that `e` and `f` are primitive idempotents, that is, idempotents in the basis.
        """        
        if e not in self.idempotents:
            raise ValueError("First argument {} is not an idempotent of {}!".format(e,self))
        if f not in self.idempotents:
            raise ValueError("Second argument {} is not an idempotent of {}!".format(f,self))

        paths = [e * x * f for x in self.basis()]
        return [p for p in paths if p != 0]
            
    # def isA1Hat(self):
    #     print("Unimplemented.")
    #     return None


    # def coeff(self, r, monomial):
    #     '''
    #     The coefficient of monomial in the expansion of r in the basis. 
    #     Monomial must be a basis element.
    #     '''
    #     i = list(self.basis()).index(monomial)
    #     return r.monomial_coefficients().get(i, 0)


    # def dualPairs(self, e, f):
    #     auxPairs = [(path, self.calabi_yau_dual(path)) for path in [x * e for x in self.basis()] if path != 0]
    #     return [(t,v*f) for (t,v) in auxPairs]
    

