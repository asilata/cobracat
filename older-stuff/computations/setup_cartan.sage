import string

load("../zigzagalgebra.sage")
load("../zigzagmodules.sage")
load("../complexes.sage")
load("../braidactions.sage")
load("../HN.sage")

q_pol.<q> = LaurentPolynomialRing(ZZ)

# Test constructor using Cartan types.
# Set up a digraph and the associated zigzag algebra, given a Cartan type.
def make_graph_and_zz_algebra(ct, k = QQ):
    """
    Given a cartan type ``ct``, return an associated digraph and a zigzag algebra.

    INPUT:
    - ``ct`` is a ``CartanType`` object, or a shorthand such as "A5" or "D7".
    
    OUTPUT:
    A tuple ``(g, Z)``, where
    - ``g`` is a labelled digraph, which is the doubled Dynkin quiver of the given Cartan type.
    - ``Z`` is the zig-zag algebra of the given Cartan type.
    """
    test_dd = CartanType(ct).dynkin_diagram()
    print(test_dd)
    test_graph_vertices = sorted(test_dd.vertices())
    test_graph_edges = []
    for (e,l) in zip(test_dd.edges(), string.ascii_lowercase):
        test_graph_edges = [(e[0], e[1], l)] + test_graph_edges
    test_graph = DiGraph([test_graph_vertices, test_graph_edges])

    return test_graph, ZigZagAlgebra(k, test_graph.path_semigroup())

def make_burau_matrices(ct) :
    """
    Given a Cartan type ``ct``, returns a dictionary ``burau_dict`` with keys the vertices of the Dynkin diagram as well as their negatives.
    For ``i`` a vertex of the Dynkin diagram, ``burau_dict[i]`` is the Burau matrix for the corresponding Artin generator, and ``burau_dict[-i]`` its inverse.

    INPUT:
    - ``ct`` a Cartan Type (or shorthand thereof)

    OUTPUT:
    A dictionary whose keys are the vertices of the Dynkin diagram as well as their negatives, and values are the corresponding Burau matrices (for the positive keys), and the inverse Burau matrices (for the negative keys).
    """
    ct = CartanType(ct)
    if not ct.is_simply_laced():
        raise NotImplementedError("{} is not a simply laced Cartan type.".format(ct))
    
    mat = matrix(CoxeterMatrix(ct)) # Computes the Coxeter matrix. 
    n = mat.nrows()
    burau_dict = {}
    #S, T = {}, {}

    for i in range(0,n) :
        loc_S = matrix(q_pol,n,n)
        loc_T = matrix(q_pol,n,n)
        for k in range(0,n) :
            # Diagonal entries
            if k == i :
                loc_S[k,k] = -q^(-2)
                loc_T[k,k] = -q^2
            else :
                loc_S[k,k] = q^0
                loc_T[k,k] = q^0

            # Other entries                
            if mat[i,k]>2 :
                loc_S[i,k] = -q^(-1)
                loc_T[i,k] = -q

                # In non simply-laced type, the formulas would be as follows
                # mult=mat[i,k]
                #loc_S[i,k] = -2*q^(-1)*cos(pi/mult)
                #loc_T[i,k] = -2*q^(-1)*cos(pi/mult)
                
            if mat[i,k]<0 : # This means that the coefficient is oo
                loc_S[i,k] = -2*q^(-1)
                loc_T[i,k] = -2*q
        burau_dict[i+1] = loc_S
        burau_dict[-(i+1)] = loc_T
    return burau_dict

# Set up the zigzag projective modules, the spherical twists, and the Burau matrices for a given Cartan Type.
def zz_setup(ct, k = QQ):
    """
    Given a cartan type ``ct``, return dictionaries of indecomposable
    projective modules, spherical twists, and inverse spherical
    twists.

    INPUT:
    - ``ct`` is a ``CartanType`` object, or a shorthand such as "A5"
      or "D7".

    OUTPUT:
    A tuple ``(p, s, t)``, where
    - ``p`` is a dictionary with keys the vertices of the Dynkin
      diagram. Given a vertex ``i``, the value ``p[i]`` is the
      associated indecomposable projective module over the zig-zag
      algebra of the cartan type ``ct``.
    - ``s`` is a dictionary with keys the vertices of the Dynkin
      diagram as well as their negatives. Given a vertex ``i``, the
      values ``s[i]`` and ``s[-i]`` are the spherical twist and
      inverse spherical twist in ``p[i]``.
    - ``burau`` is a dictionary with keys the vertices of the Dynkin
      diagram as well as their negatives. Given a vertex ``i``, the
      values ``burau[i]`` and ``burau[-i]`` are the Burau matrices
      (linearisations on the graded Grothendieck group) of ``s[i]``
      and ``s[-i]`` respectively.
    """
    graph, Z = make_graph_and_zz_algebra(ct)
    burau = make_burau_matrices(ct)

    p = {}
    for i in sorted(graph.vertices()):
        iname = "P" + str(i)
        fpi = ZigZagModule(Z, i, name = iname)
        pi = ProjectiveComplex(Z)
        pi.addObject(0, fpi)
        p[i] = pi
    
    s = {i : (lambda C, i = i: sigma(Z, i, C, minimize = True)) for i in sorted(graph.vertices())}
    s |= {-i : (lambda C, i = i: sigmaInverse(Z, i, C, minimize = True)) for i in sorted(graph.vertices())}

    return p, s, burau

# From an expression in the si's and ti's reads the terms in the Bsi's and Bti's
def braid_as_burau_list(w, burau_dict, normalize=False):
    """
    Given a word in the braid group, output a list of the Burau matrices corresponding to the generators used, from left to right.

    INPUT
    - ``w``, a word in the braid group
    - ``burau_dict``, a dictionary in the format output by ``make_burau_matrices``.
    - ``normalize``, a boolean

    OUTPUT
    A list of matrices, one per letter from left to right in the braid word.
    If ``normalize`` is True, then the braid is first converted to its left normal form.
    """
    if normalize:
        w = w.left_normal_form()[1]
        
    return [burau_dict[i] for i in w.Tietze()]

def braid_as_burau(w, burau_dict):
    """
    Given a word in the braid group, output the Burau matrix corresponding to word.

    INPUT
    - ``w``, a word in the braid group
    - ``burau_dict``, a dictionary in the format output by ``make_burau_matrices``.    

    OUTPUT
    A matrix, which is the product of the Burau matrices of the letters in the braid word.
    """
    return reduce(lambda x, y: x * y, braid_as_burau_list(w, burau_dict), identity_matrix(q_pol, max(burau_dict.keys())))

def braid_as_functor(w, s):
    """
    Given a word in the braid group, output the twist functor for the word.

    INPUT
    - ``w``, a word in the braid group
    - ``s``, a dictionary of spherical twists

    OUTPUT
    A matrix, which is the composition of the spherical twists of the letters in the braid word.
    """
    return reduce (lambda x, y: lambda t : x(y(t)), [s[i] for i in w.Tietze()], lambda t : t)

    
