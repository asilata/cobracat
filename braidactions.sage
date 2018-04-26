#load('complexes.sage')
#load('zigzagalgebra.sage')
#load('zigzagmodules.sage')

def sigmaInverse(Z, i, P):
    """
    Computes the sigma-inverse action corresponding to idempotent e_i on the projective object P. 
    The underlying zig-zag algebra is Z.
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

def sigmaInverseComplex(Z, i, C):
    e = Z.idempotents()[i-1] # Vertices are conventionally 1,2,3,... but list elements are zero-indexed :(
    Pi = ZigZagModule(Z, i, twist = 0, name="P" + str(i))
    
    # We now form a complex Q whose objects are shifts of copies of Pi 
    QObjects = {}
    mapsCtoQ = {}
    EXFs = {}
    for place in range(C.minIndex(), C.maxIndex()+1):
        QObjects[place] = {}
        mapsCtoQ[place] = {}
        EXFs[place] = {}
        for i in range(0, len(C.objects(place))):
            p = C.objects(place)[i]
            f = p.idempotent()
            EXF = filter(lambda x: x != 0, [e*x*f for x in Z.basis()])
            EXFs[place][i] = EXF
            dualPairs = [(path, Z.dualize(path)) for path in [x * e for x in Z.basis()] if path != 0]
            dualPairsf = [(t, v*f) for (t,v) in dualPairs]
            for k in range(0, len(EXF)):
                monomial = EXF[k]
                QObjects[place][(i,k)] = Pi.twistBy(p.twist() + 2 - Z.deg(monomial))
                mapsCtoQ[place][(i,k)] = sum([t * Z.coeff(Z(vf), monomial) for (t,vf) in dualPairsf])

    QMaps = {}
    for place in range(C.minIndex(), C.maxIndex()):
        QMaps[place] = {}
        for (i,k) in QObjects[place]:
            for (j, l) in QObjects[place+1]:
                r = C.maps(place).get((i,j), 0)
                targetMonomial = EXFs[place+1][j][l]
                sourceMonomial = EXFs[place][i][k]
                QMaps[place][((i,k), (j,l))] = Z.coeff(sourceMonomial * r, targetMonomial)

    # Now flatten everything
    Q = ProjectiveComplex(Z)
    keyDict = {}
    for place in range(C.minIndex(), C.maxIndex() + 1):
        keyDict[place] = {}
        counter = 0
        for (i,k) in QObjects[place].keys():
            keyDict[place][(i,k)] = counter
            counter = counter + 1
            Q.addObject(place, QObjects[place][(i,k)])

    for place in range(C.minIndex(), C.maxIndex()):
        for (i,k) in QObjects[place].keys():
            for (j,l) in QObjects[place+1].keys():
                Q.addMap(place, keyDict[place][(i,k)], keyDict[place+1][(j,l)], QMaps[place][((i,k), (j,l))])

    M = {}
    for place in range(C.minIndex(), C.maxIndex()+1):
        M[place] = {}
        for (i,k) in QObjects[place].keys():
            M[place][(i, keyDict[place][(i,k)])] = mapsCtoQ[place][(i,k)]

    print C
    print C.maps(0)
    print Q
    print Q.maps(0)
    print M
    
    return cone(C, Q, M)
