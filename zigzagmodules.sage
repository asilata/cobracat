# All modules are projective

# Left/right depends on the context. This is probably bad.

class ZigZagModule(object):
    def __init__(self, R, idempotent, twist = 0, name="P"):
        self._ring = R
        self._idempotent = idempotent
        self._twist = twist
        self._name = name


    def twistBy(self, n = 1):
        self._twist = self._twist + 1

    def _repr_(self):
        return name + "<" + str(twist) + ">"

    def _str_(self):
        return self._repr_

    def idempotent(self):
        return self._idempotent

    def tensor(self, other):
        tensors = [self._idempotent * b * other.idempotent for b in self._ring.basis()]
        answer = {}
        for t in tensors:
            answer[R.deg(t)] = answer.get(R.deg(t), 0) + 1
        return answer
