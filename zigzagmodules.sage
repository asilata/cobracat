# All modules are projective

# Left/right depends on the context. This is probably bad.

class ZigZagModule(object):
    def __init__(self, R, idempotent, twist = 0, name="P"):
        self._ring = R
        self._idempotent = idempotent
        self._twist = twist
        self._name = name


    def twistBy(self, n = 1):
        self._twist = self._twist + n

    def __repr__(self):
        return self._name + "<" + str(self._twist) + ">"

    def __str__(self):
        return self._repr_

    def idempotent(self):
        return self._idempotent

    def twist(self):
        return self._twist

    def tensor(self, other):
        tensors = [self.idempotent() * b * other.idempotent() for b in self._ring.basis()]
        nonZeroTensors = [x for x in tensors if x != 0]
        answer = {}
        for t in nonZeroTensors:
            d = R.deg(t) + self.twist() + other.twist()
            answer[d] = answer.get(d, 0) + 1
        return answer
