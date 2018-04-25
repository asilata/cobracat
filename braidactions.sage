load('complexes.sage')
load('zigzagalgebra.sage')
load('zigzagmodules.sage')

def sigma(Z, i, P):
    """
    Computes the sigma action corresponding to idempotent e on the projective object P. The underlying zig-zag algebra is Z.
    """
    e = Z.idempotents()[i-1] # Vertices are conventionally 1,2,3,... but list elements are zero-indexed :(
    f = P.idempotent()
    
    C = ProjectiveComplex(Z)
    C.addObject(0, P)
    Q = ZigZagModule(Z, i, twist=0, name="P" + str(i))

    maps = filter(lambda x: x != 0, [f*x*e for x in Z.basis()])
    if maps == []:
        return C
    else:
        for i in range(0,len(maps)):
            m = maps[i]
            C.addObject(1, Q.twistBy(Z.deg(m)))
            if Z.deg(m) == 0:
                C.addMap(0, 0, i, 1)
            else:
                C.addMap(0, 0, i, m)
        return C
    
