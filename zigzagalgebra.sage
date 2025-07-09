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
from sage.structure.element import is_Matrix

def _zz_basis(ct):
    r"""
    Return a list of basis elements of the zigzag algebra of CartanType `ct`.
    
    The standard basis of the zigzag algebra consists of the primitive
    idempotents at each vertex, arrows between two vertices, and loop
    maps at each vertex.  Each path is represented as a tuple
    `(i,j,d)`, where `i` is its source vertex, `j` is its target
    vertex, and `d` is its degree.
    
    We assume (for now) that the Dynkin diagram of the given Cartan type does not have multiple edges.
    
    INPUT:

    - `ct` -- a CartanType

    OUTPUT:

    A list of paths that form the standard basis of the zigzag algebra.

    EXAMPLES:

        sage: _zz_basis(CartanType("A2"))
        [(1, 1, 0), (2, 2, 0), (1, 2, 1), (2, 1, 1), (1, 1, 2), (2, 2, 2)]

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

def _zz_right_multiplication_table(basis, x, k=QQ):
    r"""
    Return the multiplication table for right multiplication by the basis element `x`.
    
    Generate the multiplication table for right multiplication of the elements in `basis` by a basis element `x` of the zigzag algebra.
    The element `x` must be in `basis`.
    """
    # Initialize an empty transposed matrix, which we fill in.
    mult_matrix = Matrix(k, len(basis), len(basis))
    x_source, x_target, x_deg = _zz_source(x), _zz_target(x), _zz_deg(x)

    # Case 1: x is an idempotent
    if x_deg == 0:
        for i in range(len(basis)):
            b = basis[i]
            if _zz_target(b) == x_source:
                mult_matrix[i,i] = 1

    # Case 2: x is an edge (degree-one element) 
    elif x_deg == 1:
        for i in range(len(basis)):
            b = basis[i]
            if _zz_deg(b) == 0 and _zz_target(b) == x_source:
                # Multiplication with idempotent at source of x is x.
                j = basis.index(x)
                mult_matrix[i,j] = 1
            elif (_zz_deg(b) == 1 and
                  _zz_target(b) == x_source and 
                  _zz_source(b) == x_target):
                # Multiplication with the opposite edge produces the loop at the target of x.
                # CHECK SIGN CONVENTION HERE.
                j = basis.index((x_target, x_target,2))
                if _zz_source(b) < _zz_target(b):
                    # b is a "forwards" arrow
                    mult_matrix[i,j] = 1
                else:
                    # b is a "backwards" arrow
                    mult_matrix[i,j] = -1                    
    # Case 3: x is a loop
    elif x_deg == 2:
        for i in range(len(basis)):
            b = basis[i]
            if _zz_deg(b) > 0 or (_zz_target(b) != x_source):
                # Multiplication with positive degree elements as well as
                # idempotents that don't match the source of x is zero.
                pass
            else:
                # Multiplication with the idempotent with the same source and target yields x.
                j = basis.index(x)
                mult_matrix[i,j] = 1
    else:
        raise ValueError("{} is not a basis element of the zigzag algebra!".format(x))
    return mult_matrix

def _ct_has_multiple_edges(ct):
    return not all([all([y >= 0 or y == -1 for y in r]) for r in ct.cartan_matrix()])

class ZigZagAlgebra(FiniteDimensionalAlgebra):
    r"""
    The zig-zag algebra over a field associated to a Cartan type.
    """
    Element = ZigZagAlgebraElement
    
    def __init__(self, ct, k=QQ, index_set=None):
        r"""
        Create a zig-zag algebra.

        INPUT:

        - `ct` -- A matrix, a `CartanMatrix`, a `CartanType`, or a shorthand such as "A5" or "D4". A shorthand input as a list (e.g. ['A', 5, 1]) is not supported.
        - `k` -- Field, default `QQ`
        - `index_set` -- a tuple to reindex the vertices of the Dynkin diagram, default `None`.

        OUTPUT:

        The zig-zag algebra of type `ct` over the field `k`.  The zig-zag algebra is
        `Associative`, `Graded`, `FiniteDimensional`, `Algebras(k)`, and `WithBasis`.
        
        """

        if is_Matrix(ct):
            if index_set is None:
                index_set = range(1,ct.nrows()+1)
            ct = CartanMatrix(ct, index_set)
            
        # if isinstance(ct, sage.combinat.root_system.cartan_matrix.CartanMatrix):
        #     if index_set is None:
        #         index_set = range(1,ct.cartan_matrix().nrows()+1)
        #     if len(index_set) != ct.cartan_matrix().nrows():
        #         raise ValueError("Size of index set {} does not match size of Cartan matrix!".format(index_set))
        #     ct = DynkinDiagram(ct, index_set=index_set)
            
        elif not isinstance(ct, sage.combinat.root_system.cartan_type.CartanType_abstract):
            ct = CartanType(ct)

        self.cartan_type = ct
        self.index_set = ct.index_set()
            
        self.cartan_matrix = self.cartan_type.cartan_matrix()
            
        if not self.cartan_matrix.is_symmetric():
            raise ValueError("Cartan matrix not symmetric: {}".format(self.cartan_matrix))
        
        if _ct_has_multiple_edges(self.cartan_type):
            raise NotImplementedError("ZigZagAlgebra not implemented for Dynkin diagrams with multiple edges.")
            
        self._base_ring = k

        # An internal representation for the basis of `self`, as specified in the helper methods `_zz_basis()`defined previously.
        self._internal_basis = _zz_basis(self.cartan_type)

        # A list of matrices, where the ith matrix is the matrix of right multiplication
        # by the ith basis element
        table = [_zz_right_multiplication_table(self._internal_basis, b, k) for b in self._internal_basis]
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

    def calabi_yau_dual(self, b):
        r"""
        Return the unique Calabi--Yau dual basis element of `b`.

        INPUT:

        - `b`  -- a basis element of `self`.

        OUTPUT:

        - The unique basis element `a` of `self` such that `a * b` and `b * a` are both plus or minus of the loops at the corresponding vertices.
          Raise `ValueError` if `b` is not in `self.basis()`.
        """
        if b not in self.basis():
            raise ValueError("{} is not a basis element of {}.".format(b,Z))
        s,t = self.idempotent_by_vertex(self.source(b)), self.idempotent_by_vertex(self.target(b))
        for c in self.basis():
            bc = b * t * c * s
            if bc in self.loops or -bc in self.loops:
                return c

    @property
    def cy_duals(self):
        output = []
        for b in self.basis():
            output.append((b, self.calabi_yau_dual(b)))
        return output        

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
        Return the source vertex of the given basis element.

        INPUT:

        - `b` -- A basis element of `self`.

        OUTPUT:

        The source vertex of the arrow represented by `b`.  Raise `ValueError` if `b` is not an element of `self.basis()`.
        """
        if b not in self.basis():
            raise ValueError("{} is not a basis element of {}!".format(b,self))

        return _zz_source(self._to_internal_repr(b))

    def target(self,b):
        """
        Return the target vertex of the given basis element.

        INPUT:

        - `b` -- A basis element of `self`.

        OUTPUT:

        The target vertex of the arrow represented by `b`.  Raise `ValueError` if `b` is not an element of `self.basis()`.
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

    @cached_method
    def idempotent_by_vertex(self, v):
        """
        Returns the idempotent of `self` corresponding to vertex `v`.
        """
        if v not in self.index_set:
            raise ValueError("{} is not a vertex of {}!".format(v,self))
        # Return the first (and only) idempotent that has v as a source vertex.
        return [e for e in self.idempotents if self.source(e) == v][0]

    @cached_method
    def coeff(self, r, monomial):
        '''
        The coefficient of monomial in the expansion of r in the basis. 
        Monomial must be a basis element.
        '''
        i = list(self.basis()).index(monomial)
        return r.monomial_coefficients().get(i, 0)

    @cached_method
    def paths_from_to(self,e,f):
        r"""
        Returns a list of basis elements that represent paths from the source=target of `e` to the source=target of `f`.

        INPUT:

        - `e` -- a primitive idempotent, that is, an idempotent in `self.basis()`
        - `f` -- a primitive idempotent, that is, an idempotent in `self.basis()`

        OUTPUT:

        - A list of basis elements each of which represents a path from the source=target of `e` to the source=target of `f`.
          Raise `ValueError` if `e` or `f` is not an idempotent in `self.basis()`.
        """        
        if e not in self.idempotents:
            raise ValueError("First argument {} is not an idempotent of {}!".format(e,self))
        if f not in self.idempotents:
            raise ValueError("Second argument {} is not an idempotent of {}!".format(f,self))

        paths = [e * x * f for x in self.basis()]
        return [p for p in paths if p != 0]

    

