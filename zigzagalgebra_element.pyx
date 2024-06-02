from sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra_element cimport FiniteDimensionalAlgebraElement

cdef class ZigZagAlgebraElement(FiniteDimensionalAlgebraElement):
    def __init__(self, A, elt=None, check=True):
        FiniteDimensionalAlgebraElement.__init__(self, A = A, elt = elt, check = check)
        

