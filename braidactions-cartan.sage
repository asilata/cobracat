def sigma(Z, i, C):
    '''
    The spherical twist corresponding to the i-th projective module of
    the ZigZagAlgebra Z.

    More explicitly, let Pi = Z * Z.idempotent_by_vertex(i). Consider the map
    K of bimodules Z -> Pi tensor_k iP. Then sigma(Z, i, C) is the
    cone over the map of complexes obtained by tensoring K by C.

    sigma(Z, i, C) = Cone(P_i tensor_k iP tensor C -> C)[1]

    '''
    e = Z.idempotent_by_vertex(i)
    # if Z.isA1Hat():
    #     Pi = A1HatModule(Z, i, twist = 0, name="P" + str(i))
    # else:
    #     Pi = ZigZagModule(Z, i, twist = 0, name="P" + str(i))
    Pi = ProjectiveZigZagModule(Z, i, graded_degree = 0, name_prefix = "P")

    # We now form a complex Q whose objects are shifts of copies of Pi 
    QObjects = {}
    mapsQtoC = {}
    EXFs = {}
    for place in range(C.min_index, C.max_index+1):
        QObjects[place] = {}
        mapsQtoC[place] = {}
        EXFs[place] = {}
        for i in range(0, len(C.objects[place])):
            p = C.objects[place][i]
            f = p.idempotent
            EXF = Z.paths_from_to(e,f)
            EXFs[place][i] = EXF
            for k in range(0, len(EXF)):
                monomial = EXF[k]
                QObjects[place][(i,k)] = Pi.graded_shift_by(p.graded_degree - Z.deg(monomial))
                mapsQtoC[place][(i,k)] = monomial

    QMaps = {}
    for place in range(C.min_index, C.max_index):
        QMaps[place] = {}
        for (i,k) in QObjects[place]:
            for (j, l) in QObjects[place+1]:
                r = C.maps[place].get((i,j), 0)
                targetMonomial = EXFs[place+1][j][l]
                sourceMonomial = EXFs[place][i][k]
                QMaps[place][((i,k), (j,l))] = Z.coeff(sourceMonomial * r, targetMonomial)

    # Now flatten everything
    Q = ProjectiveComplex(Z)
    keyDict = {}
    for place in range(C.min_index, C.max_index + 1):
        keyDict[place] = {}
        counter = 0
        for (i,k) in QObjects[place].keys():
            keyDict[place][(i,k)] = counter
            counter = counter + 1
            Q.add_object_at(place, QObjects[place][(i,k)])

    for place in range(C.min_index, C.max_index):
        for (i,k) in QObjects[place].keys():
            for (j,l) in QObjects[place+1].keys():
                Q.add_map_at(place, keyDict[place][(i,k)], keyDict[place+1][(j,l)], QMaps[place][((i,k), (j,l))])

    M = {}
    for place in range(C.min_index, C.max_index+1):
        M[place] = {}
        for (i,k) in QObjects[place].keys():
            M[place][(keyDict[place][(i,k)], i)] = mapsQtoC[place][(i,k)]

    return cone(Q, C, M).homological_shift_by(1)
    

def sigmaInverse(Z, i, C):
    '''The inverse spherical twist corresponding to the i-th projective
    module of the ZigZagAlgebra Z.

    More explicitly, let Pi = Z * Z.idempotents()[i]. Consider the map
    K of bimodules Z -> Pi tensor_k iP. Then sigmaInverse(Z, i, C) is
    the cone over the map of complexes obtained by tensoring K by C.

    sigmaInverse(Z, i, C) = Cone(C -> P_i tensor_k iP tensor C)

    '''
    e = Z.idempotent_by_vertex(i)
    # if Z.isA1Hat():
    #     Pi = A1HatModule(Z, i, twist = 0, name="P" + str(i))
    # else:
    #     Pi = ZigZagModule(Z, i, twist = 0, name="P" + str(i))
    Pi = ProjectiveZigZagModule(Z, i, graded_degree = 0, name_prefix = "P")

    # We now form a complex Q whose objects are shifts of copies of Pi 
    QObjects = {}
    mapsCtoQ = {}
    EXFs = {}
    for place in range(C.min_index, C.max_index+1):
        QObjects[place] = {}
        mapsCtoQ[place] = {}
        EXFs[place] = {}
        for i in range(0, len(C.objects[place])):
            p = C.objects[place][i]
            f = p.idempotent
            EXF = Z.paths_from_to(e,f)
            EXFs[place][i] = EXF
            for k in range(0, len(EXF)):
                monomial = EXF[k]
                QObjects[place][(i,k)] = Pi.graded_shift_by(p.graded_degree + 2 - Z.deg(monomial))
                mapsCtoQ[place][(i,k)] = sum([t * Z.coeff(Z(vf), monomial) for (t,vf) in Z.dualPairs(e,f)])

    QMaps = {}
    for place in range(C.min_index, C.max_index):
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
    for place in range(C.min_index, C.max_index + 1):
        keyDict[place] = {}
        counter = 0
        for (i,k) in QObjects[place].keys():
            keyDict[place][(i,k)] = counter
            counter = counter + 1
            Q.add_object_at(place, QObjects[place][(i,k)])

    for place in range(C.min_index, C.max_index):
        for (i,k) in QObjects[place].keys():
            for (j,l) in QObjects[place+1].keys():
                Q.add_map_at(place, keyDict[place][(i,k)], keyDict[place+1][(j,l)], QMaps[place][((i,k), (j,l))])

    M = {}
    for place in range(C.min_index, C.max_index+1):
        M[place] = {}
        for (i,k) in QObjects[place].keys():
            M[place][(i, keyDict[place][(i,k)])] = mapsCtoQ[place][(i,k)]

    return cone(C, Q, M)

from functools import reduce
def composeAll(list_of_functions):
    return reduce (lambda x,y : lambda t : x(y(t)), list_of_functions, lambda x: x)
