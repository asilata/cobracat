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

    Element = ZigZagModuleElement
    
    def __init__(self, Z, v, graded_degree = 0, name_prefix='P'):
        '''
        Initialize `self`.

        `Z` is a ZigZagAlgebra, and v is a vertex of Z. The `graded_degree` is the internal graded degree, which defaults to zero.
        `name_prefix` is the prefix for the names, (default 'P').
        The output is the projective module P = Z*ev, where ev is the idempotent corresponding to vertex v in Z.
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
        
