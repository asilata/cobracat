from functools import reduce

def sigma_new(C, i, Z, minimize=False):
    r"""
    The spherical twist corresponding to the ith projective module of the zig-zag algebra Z.

    More explicitly, let H be the complex Hom(Pi, C).
    We have a universal evaluation map from Hom(Pi, C) tensor Pi to C.
    We construct this chain map, and then return its cone.
    """

    e = Z.idempotent_by_vertex(i)
    Pi = ProjectiveZigZagModule(Z, i, graded_degree=0, name_prefix="P")
    Pi_complex = ProjectiveComplex(Z)
    Pi_complex.add_object_at(0, Pi)
    
    H, chain_map_data = Pi_complex.hom(C, with_homs=True)
    # Assume that H = Hom(Pi_complex, C) is already minimized.

    # Construct the complex H tensor Pi.
    # All objects of this complex are copies of Pi with homological and degree shifts.
    # This complex has no differentials internally, because H has no differentials.
    A = ProjectiveComplex(algebra=Z)
    for j in H.objects:
        for k_d in H.objects[j]:
            A.add_object_at(j, Pi.graded_shift_by(k_d.graded_degree))

    # Construct the chain map from A to C.
    # This is constructed using `chain_map_data`.
    for x in chain_map_data:
        for y in chain_map_data[x]:
            assert set(chain_map_data[x][y].keys()) == set([0])

    M = {}
    for homological_index in H.objects:
        if homological_index not in M:
            M[homological_index] = {}
            
        for a in range(len(H.objects[homological_index])):
            chain_map_a = chain_map_data[homological_index][a]

            # chain_map_a is a chain map from Pi_complex, which is
            # concentrated in degree 0, to C. The only possible key it can have is 0.
            assert 0 in chain_map_a and len(chain_map_a) == 1

            # Furthermore, chain_map_a[0] only contains entries of the form
            # (0, b), signifying a chain map to the bth object of C.
            assert all([c == 0 for (c,_) in chain_map_a[0]])

            # Now for each (0,b) in chain_map_a[0], add the map chain_map_a[0][(0,b)] to A.
            for (c,b) in chain_map_a[0]:
                M[homological_index][(a,b)] = chain_map_a[0][(c,b)]
                
    answer = cone(A, C, M).homological_shift_by(1)
    if minimize:
        return answer.minimize_using_matrix()
    else:
        return answer

def sigma(Z, i, C, minimize=False):
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

    # TODO the following code assumes that C is "clean"; that is, that
    # C.cleanup() does not change it. Make more robust to remove this assumption.
    
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

    c = cone(Q, C, M).homological_shift_by(1)
    if minimize:
        c = c.minimize_using_matrix()
    return c

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

    c = cone(C,Q,M)
    if minimize:
        c = c.minimize_using_matrix()
    return c

def composeAll(list_of_functions):
    return reduce (lambda x,y : lambda t : x(y(t)), list_of_functions, lambda x: x)
