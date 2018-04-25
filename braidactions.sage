load('complexes.sage')
load('zigzagalgebra.sage')
load('zigzagmodules.sage')

def sigmaInverse(Z, i, P):
    """
    Computes the sigma-inverse action corresponding to idempotent e on the projective object P. The underlying zig-zag algebra is Z.
    """
    e = Z.idempotents()[i-1] # Vertices are conventionally 1,2,3,... but list elements are zero-indexed :(
    f = P.idempotent() # So P = Pf.
    
    C = ProjectiveComplex(Z)
    C.addObject(0, P) # No matter where e is with respect to f, we always first add a copy of P in the 0th place.
    Q = ZigZagModule(Z, i, twist=0, name="P" + str(i)) # This is Pe.

    maps = filter(lambda x: x != 0, [f*x*e for x in Z.basis()]) # Figure out how many maps there are from Pf to Pe.
    if maps == []:
        return C # If no maps, this means e is far away from f. Return C as-is.
    else:
        for i in range(0,len(maps)):
            m = maps[i]
            C.addObject(1, Q.twistBy(Z.deg(m))) # For each arrow from f to e, add a Pe with the right degree twist.
            C.addMap(0, 0, i, m) # Add the appropriate map Pf to Pe<twist>.
        return C
    
