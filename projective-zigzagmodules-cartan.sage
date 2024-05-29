r"""
Projective (left) modules over zig-zag algebras

We implement indecomposable left projective modules over zig-zag algebras. Given a zig-zag algebra Z and a vertex v of the quiver of Z, we have the primitive idempotent e, which is the length-zero path at v. We construct the left module Ze. This module has a natural grading inherited from Z.

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

from sage.modules.module import Module
from sage.structure.element import ModuleElement

class ZigZagModuleElement(ModuleElement):
    def __init__(self, parent, x):
        self.x = x
        ModuleElement.__init__(self, parent = parent)

class ProjectiveZigZagModule(Module):
    '''
    The indecomposable projective module over a given zigzag algebra, corresponding to a given vertex.
    '''

    # We intern the zig-zag modules.
    # All created instances are stored here.
    _instances = {}

    # First check if an instance with the given parameters has been created.
    # If yes, return that.  Otherwise, create a new one.
    def __new__(cls, Z, v, graded_degree = 0, name_prefix='P'):
        r'''
        TESTS:

        sage: load("zigzagalgebra-cartan.sage")
        sage: Z = ZigZagAlgebra("A2", QQ)
        sage: P0 = ProjectiveZigZagModule(Z, Z.vertices[0])
        sage: P0p = ProjectiveZigZagModule(Z, Z.vertices[0])
        sage: P0 == P0p
        True
        sage: P1 = ProjectiveZigZagModule(Z, Z.vertices[1])
        sage: P0 == P1
        False
        '''
        if (Z,v,graded_degree,name_prefix) in cls._instances:
            return cls._instances[(Z,v,graded_degree,name_prefix)]
        else:
            instance = super().__new__(cls)
            cls._instances[(Z,v,graded_degree,name_prefix)] = instance
            return instance
        
    Element = ZigZagModuleElement

    
    def __init__(self, Z, v, graded_degree = 0, name_prefix='P'):
        '''
        Initialize `self`.

        INPUT: 

        - `Z` -- a `ZigZagAlgebra`
        - `v` --- a vertex of `Z`
        - `graded_degree` -- integer, the internal graded degree (default 0)
        - `name_prefix` -- string, a prefix for the names (default 'P')

        OUTPUT:

        - The projective module `P` = `Z * ev`, grade-shifted by `graded_degree`.
          Here `ev` is the idempotent corresponding to the vertex `v` in `Z`.
        '''
        self.algebra = Z
        self.vertex = v
        self.idempotent = Z.idempotent_by_vertex(v)
        self.graded_degree = graded_degree
        self._name_prefix = name_prefix
        Module.__init__(self, Z, category=LeftModules(Z))

    def __repr__(self):
        return self._name_prefix + str(self.vertex) + "<" + str(self.graded_degree) + ">"

    def __str__(self):
        return self.__repr__()
        
    def copy(self):
        return ZigZagModule(self.algebra, self.vertex, self.graded_degree, self._name_prefix)
        
    def graded_shift_by(self, n = 1):
        """
        Return a copy of `self` that is internally degree shifted by n.        
        """
        return ZigZagModule(self.algebra, self.vertex, self.graded_degree + n, self._name_prefix)

    #@cached_method
    def hom(self, Q):
        """
        Return a list of elements of the zigzag algebra that send `self` to Q on right multiplication.
        """
        e = self.idempotent
        f = Q.idempotent
        return [e*x*f for x in self.algebra.basis() if e*x*f != 0]

    # # Here the module is supposed to be the left R-module Re, but the map is right multiplication by r.
    # @cached_method
    def is_annihilated_by(self, r):
        r"""
        Returns True if right multiplication by `r` is the zero map on `self`.
        """        
        return (self.idempotent * r == 0)

    @cached_method
    def is_invertible(self, r):
        r"""
        Returns True if right multiplication by `r` is an invertible map on `self`.
        """        
        ir = self.idempotent * r
        if ir == 0:
            return False
        non_zero_coeffs = [(x,y) for (x,y) in self.idempotent.monomial_coefficients().items() if y != 0]
        m,d = non_zero_coeffs[0]
        c = ir.monomial_coefficients().get(m, 0)
        multiple = c/d
        return (ir == multiple*self.idempotent)

    @cached_method
    def invert(self, r):
        r"""
        Returns an element s in basering such that right multiplication by s is the inverse of right multiplication by r.
        Assumes that right multiplication by r is invertible.
        """        
        if not self.is_invertible(r):
            raise TypeError("Not invertible")
        else:
            ir = self.idempotent * r
            non_zero_coeffs = [(x,y) for (x,y) in self.idempotent.monomial_coefficients().items() if y != 0]
            m,d = non_zero_coeffs[0]
            c = ir.monomial_coefficients().get(m, 0)
            multiple = c/d
            return d/c * self.idempotent
        
