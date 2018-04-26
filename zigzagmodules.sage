# All modules are projective

# Left/right depends on the context. This is probably bad.

class ZigZagModule(object):
    def __init__(self, R, i, twist = 0, name="P"):
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

    def tensor(self, other):
        tensors = [self.idempotent() * b * other.idempotent() for b in self._ring.basis()]
        nonZeroTensors = [x for x in tensors if x != 0]
        answer = {}
        for t in nonZeroTensors:
            d = R.deg(t) + self.twist() + other.twist()
            answer[d] = answer.get(d, 0) + 1
        return answer

    # Here the module is supposed to be the left R-module Ri, but the map is right multiplication by r.
    def is_zero(self, r):
        return (self._idempotent * r == 0)

    def is_invertible(self, r):
        ir = self._idempotent * r
        if ir == 0:
            return False

        c = ir.coefficients(sort=True)[0]
        d = self._idempotent.coefficients(sort=True)[0]
        multiple = c/d

        return (ir == multiple*self._idempotent)

    def invert(self, r):
        if not self.is_invertible(r):
            raise TypeError("Not invertible")
        else:
            ir = self._idempotent * r
            c = ir.coefficients(sort=True)[0]
            d = self._idempotent.coefficients(sort=True)[0]
            multiple = c/d
            return d/c * self._idempotent
        
