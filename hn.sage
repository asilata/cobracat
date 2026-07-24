def find_top(obj, stab):
    """
    Compute the top HN piece of the given object with respect to the stability condition being considered.

    INPUT:

    - `obj` : a ProjectiveComplex
    - `stab`: a sequence of objects ordered by increasing phase.

    OUTPUT:

    A tuple (x, h) where x is a stable object of the highest phase that has a non-trivial morphism to obj.
    The morphism h is one such morphism.
    """
    highest_t_deg_found = -Infinity
    highest_object_found = None
    highest_hom_found = None

    for stable in stab:
        rawHoms = stable.hom(obj, unreduced=True, collapsed_grading=True)
        for i in range(rawHoms.max_index, rawHoms.minindex()-1,-1):
            betti = len(rawHoms[i]) - matrix(rawHoms.maps.get(i,{})).rank() - matrix(rawHoms.maps.get(i-1,{})).rank()
            if betti > 0 and i >= highest_t_deg_found:
                highest_t_deg_found = i
                highest_obj_found = stable
    
    top = highest_obj_found
    (h, hom_data)= top.hom(obj, with_homs=True)

    highest_t_deg_found = -Infinity
    highest_object_found = None
    highest_hom_found = None
    for i in range(h.min_index, h.max_index+1):
        for j in range(0, len(h[i])):
            if h[i][j].graded_degree - i >= highest_t_deg_found:
                highest_t_deg_found = h[i][j].graded_degree - i
                highest_object_found = stable.graded_shift(h[i][j].graded_degree).homological_shift(-i)
                highest_hom_found = {x+i : hom_data[i][j][x] for x in hom_data[i][j]}
                    
    return (highest_object_found, highest_hom_found)

def hn_factors(obj, stab):
    """
    Compute the HN filtration of object obj with respect to the given (generic) stability condition.

    INPUT:

    - `obj` : a ProjectiveComplex.
    - `stab`: a sequence of objects ordered by increasing phase

    OUTPUT:

    A list of stable HN (or refined HN) factors of obj, in non-increasing phase order.
    """
    hn_factors_list = []
    
    while not obj.is_zero():
        (top_piece, top_hom) = find_top(obj, stab)
        hn_factors_list.append(top_piece)
        obj = cone(top_piece, obj, top_hom).homological_shift(1)
        obj.minimize()

    return hn_factors_list


