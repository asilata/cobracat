from functools import reduce

def sigma(C, i, Z, minimize=False):
    r"""
    The spherical twist corresponding to the ith projective module of the zig-zag algebra Z.

    More explicitly, let H be the complex Hom(Pi, C).
    We have a universal evaluation map from Hom(Pi, C) tensor Pi to C.
    We construct this chain map, and then return its cone, homologically shifted by [1].
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
    M = {}
    for homological_index in H.objects:
        if homological_index not in M:
            M[homological_index] = {}
            
        for a in range(len(H.objects[homological_index])):
            chain_map_a = chain_map_data[homological_index][a]

            # chain_map_a is a chain map from Pi_complex, which is
            # concentrated in degree 0, to C. The only possible key it can have is 0.
            # This assertion assumes that chain_map_a does not have spurious zero entries.
            assert 0 in chain_map_a and len(chain_map_a) == 1

            # Furthermore, chain_map_a[0] only contains entries of the form
            # (0, b), signifying a chain map to the bth object of C.
            assert all([c == 0 for (c,_) in chain_map_a[0]])

            # Now for each (0,b) in chain_map_a[0], add the map chain_map_a[0][(0,b)] to M.
            for (c,b) in chain_map_a[0]:
                M[homological_index][(a,b)] = chain_map_a[0][(c,b)]
                
    answer = cone(A, C, M).homological_shift_by(1)
    if minimize:
        return answer.minimize_using_matrix()
    else:
        return answer

def sigma_inverse(C, i, Z, minimize=False):
    r"""
    The inverse spherical twist corresponding to the ith projective module of the zig-zag algebra Z.

    More explicitly, let H be the dual of the complex Hom(C, Pi).
    We have a universal co-evaluation map from C to H tensor Pi.
    We construct this chain map, and then return its cone.
    """
    e = Z.idempotent_by_vertex(i)
    Pi = ProjectiveZigZagModule(Z, i, graded_degree=0, name_prefix="P")
    Pi_complex = ProjectiveComplex(Z)
    Pi_complex.add_object_at(0, Pi)
    
    Hdual, chain_map_data = C.hom(Pi_complex, with_homs=True)
    # Assume that Hdual is already minimized.

    # Construct the complex H tensor Pi
    # All objects of this complex are copies of Pi with homological and degree shifts.
    # This complex has no differentials internally, because Hdual has no differentials.
    A = ProjectiveComplex(algebra=Z)
    for j in Hdual.objects:
        for a in range(len(Hdual.objects[j])):
            k_d = Hdual.objects[j][a]
            A.add_object_at(-j, Pi.graded_shift_by(-k_d.graded_degree))

    # Construct the chain map from C to A
    # This is constructed using `chain_map_data`
    M = {}
    for homological_index in Hdual.objects:
        if -homological_index not in M:
            M[-homological_index] = {}
        for a in range(len(Hdual.objects[homological_index])):
            chain_map_a = chain_map_data[homological_index][a]

            # Let d be the internal degree of this map.
            # The source of chain_map_a is C
            # The target of chain_map_a is Pi<d>[homological_index]
            # So, if we let k<-d>, living in homological degree
            # homological_index, be the summand corresponding to this
            # map, then the co-evaluation map should be
            # chain_map_a: C -> k<d> tensor Pi
            # the target living in homological degree -homological_index.
            
            # This assertion assumes that chain_map_a does not have spurious zero entries.
            assert -homological_index in chain_map_a and len(chain_map_a) == 1
            
            # Furthermore, chain_map_a[-homological_index] only
            # contains entries of the form (c, 0) signifying a map
            # from the cth object of C to the 0th and only object of
            # Pi[homological_index].
            assert all([b == 0 for (_,b) in chain_map_a[-homological_index]])

            # For each (c,0) in chain_map_a[-homological_index], add the map
            # chain_map_a[-homological_index][(c,0)] to M
            for (c,b) in chain_map_a[-homological_index]:
                M[-homological_index][(c,a)] = chain_map_a[-homological_index][(c,b)]
                
    answer = cone(C, A, M)
    if minimize:
        return answer.minimize_using_matrix()
    else:
        return answer

def composeAll(list_of_functions):
    return reduce (lambda x,y : lambda t : x(y(t)), list_of_functions, lambda x: x)
