class ZigZagModule(object):
'''
Projective (left) modules over the Zigzag algebra
'''
    def __init__(self, R, i, twist = 0, name="P"):
        '''
        The projective module P = R*e, where e is the i-th idempotent in R.idempotents()
        '''
        self._ring = R
        self._i = i
        self._idempotent = R.idempotents()[i-i]  # Vertices are conventionally 1,2,3,... but list elements are zero-indexed :(
        self._twist = twist
        self._name = name

    def copy(self):
        return ZigZagModule(self._ring, self._i, self._twist, self._name)
        
    def twistBy(self, n = 1):
        return ZigZagModule(self._ring, self._i, self._twist + n, self._name)

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

    # Here the module is supposed to be the left R-module Re, but the map is right multiplication by r.
    def is_zero(self, r):
        return (self._idempotent * r == 0)

    def is_invertible(self, r):
        ir = self._idempotent * r
        if ir == 0:
            return False

        m,d = self._idempotent.monomial_coefficients().items()[0]
        c = ir.monomial_coefficients().get(m, 0)
        multiple = c/d

        return (ir == multiple*self._idempotent)

    def invert(self, r):
        if not self.is_invertible(r):
            raise TypeError("Not invertible")
        else:
            ir = self._idempotent * r
            m,d = self._idempotent.monomial_coefficients().items()[0]
            c = ir.monomial_coefficients().get(m, 0)
            multiple = c/d
            return d/c * self._idempotent
        
