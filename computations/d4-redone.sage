# load("../complexes.sage")
# load("../zigzagalgebra.sage")
# load("../zigzagmodules.sage")
# load("../braidactions.sage")

import itertools
import networkx as nx

load("setup_cartan.sage")

# Set the Cartan type here.
ct = CartanType("D4")

# p = dictionary of projectives indexed by the vertices of the Dynkin diagram
# s = dictionary of positive and negative spherical twists, indexed by the vertices of the Dynkin diagram and their negatives.
# burau = dictionary of Burau matrices, indexed by the vertices of the Dynkin diagram and their negatives.
p, s, burau = zz_setup(ct)

# Weyl group of the chosen Cartan type
W = WeylGroup(ct)

def all_signed_combinations(w):
    """
    Return a list of all possible signed lifts of a given word w in the Weyl group.

    INPUT:
    - w, a word in the Weyl group, for example [1,2]

    OUTPUT:
    - a list of all possible signed versions of the given word. In the previous example that would be [[1,2], [-1,2], [1,-2], [-1,-2]].
    """
    current_signed_lifts = [[]]
    current_word = w

    while(len(current_word) > 0):
        new_signed_lifts = []
        current_letter = current_word[0]
        current_word = current_word[1:]
        for v in current_signed_lifts:
            new_signed_lifts.append(v + [current_letter])
            new_signed_lifts.append(v + [-current_letter])
        current_signed_lifts = new_signed_lifts
    return current_signed_lifts

# A sorted list of all the reduced words in the Weyl group, written as lists in the vertices of the Dynkin diagram.
reduced_words = sorted([x.reduced_word() for x in W], key=len)

# Collect all signed combinations of all reduced words in the Weyl group.
reduced_braid_words = sum([all_signed_combinations(w) for w in reduced_words], [])

def in_heart(obj):
    """
    Return True if and only if the given object is in the standard heart; that is, if it is a linear complex.
    """
    for i in range(obj.minIndex(), obj.maxIndex()+1):
        obj_i = obj.objects(i)
        for x in obj_i:
            if x.twist() != i:
                return False
    return True
    
def equal_upto_shift(x,y):
    """
    Check if objects x and y are equal up to shift.
    """
    xmin, xmax = x.minIndex(), x.maxIndex()
    ymin, ymax = y.minIndex(), y.maxIndex()
    if xmax - xmin != ymax - ymin:
        return False

    for i in range(0, xmax - xmin + 1):
        xi, yi = x.objects(i + xmin), y.objects(i + ymin)
        if len(xi) !=  len(yi):
            return False
        if sorted([a.name() for a in xi]) != sorted([a.name() for a in yi]):
            return False
    return True

def generate_elements_in_heart():
    """
    Return a list of all the spherical objects in the standard heart of the 2-CY category associated to the chosen Cartan type.
    This list checks if the generated object at any stage is equal (up to shift) to any previously generated object, and only adds it to the outputs if not already present.
    Thus the returned list should contain no duplicates.
    """
    # Start with the list of reduced braid words; this will be pruned as we go.
    pruned_braid_words = reduced_braid_words
    
    long_word_length = max([len(x) for x in reduced_words])    
    outputs = []

    # Go through words by length, applying the word to p[1] and checking if it is in the heart.
    # If it is in the heart, then add it to outputs.
    # Otherwise, prune the remaining braid words by removing any word that has the current word as a rightmost substring.
    for l in range(0,long_word_length + 1):
        print("Applying words of length {}.".format(l))
        lifts_of_length_l = [x for x in pruned_braid_words if len(x) == l]
        pruned_braid_words = [x for x in pruned_braid_words if len(x) > l]
        
        for b in lifts_of_length_l:
            b_of_p1 = composeAll([s[i] for i in b])(p[1])
            if in_heart(b_of_p1):
                exist_equal_elements = False
                for x in outputs:
                    if equal_upto_shift(x,b_of_p1):
                        exist_equal_elements = True
                        break

                if not exist_equal_elements:
                    outputs = outputs + [b_of_p1]
            else:
                pruned_braid_words = [x for x in pruned_braid_words if x[-len(b):] != b]
    return outputs

def exists_hom_one(A,B):
    """
    Return true if there is at least one degree-one hom from object A to object B.
    """
    h = hom(A,B).qPolynomial()
    hom_ones = [c for c in h.coefficients() if c[1] == -1]
    if len(hom_ones) >= 1:
        return True
    else:
        return False

def find_no_hom_one_cliques(candidates):
    """
    Given a list of candidate objects, return maximum-size cliques of objects that pairwise have no degree one homs between them.
    """
    no_hom_one_pairs = [(candidates[i], candidates[j]) for i in range(0, len(candidates)) for j in range(0,i) if not exists_hom_one(candidates[i], candidates[j])]
    G = nx.Graph()
    G.add_edges_from(no_hom_one_pairs)
    g_cliques = list(nx.find_cliques(G))
    m = max([len(x) for x in g_cliques])
    return [x for x in g_cliques if len(x) == m]

def no_hom_one_graph(candidates):
    no_hom_one_pairs = [(candidates[i], candidates[j]) for i in range(0, len(candidates)) for j in range(0,i) if not exists_hom_one(candidates[i], candidates[j])]
    G = nx.Graph()
    G.add_edges_from(no_hom_one_pairs)
    G = Graph(G)
    return G

def find_heart_and_ppts():
    """
    Return a list of all spherical objects in the heart, as well as all the ppts in the heart, defined as maximum cliques of objects that have pairwise no degree-one homs.
    """
    heart = generate_elements_in_heart()
    ppts = [frozenset(x) for x in find_no_hom_one_cliques(heart)]
    return heart, ppts

def flip_graph(ppts):
    max_intersection_pairs = [(ppts[i], ppts[j]) for i in range(0, len(ppts)) for j in range(0,i) if len(ppts[i] & ppts[j]) == len(ppts[i]) - 1]
    G = nx.Graph()
    G.add_edges_from(max_intersection_pairs)
    return G

# Graveyard beyond this point.

def find_external_edges(candidates):
    externals = []
    for x in candidates:
        hom_ones = [y for y in candidates if exists_hom_one(x,y)]
        if len(hom_ones) > 0:
            continue
        else:
            externals = externals + [x]
    return externals

def find_flips(p, ppts):
    flip_dictionary = {}
    for e in p:
        flips = [q for q in ppts if q != p and q - set([e]) == p - set([e])]
        flip_dictionary[e] = flips
    return flip_dictionary

def flip_edge(e, ppt, ppts):
    set_int = set(ppt) - set([e])
    return [q for q in ppts if set(q) & set(ppt) == set_int]

def flip_all(ppts):
    unflippable = {}
    for p in ppts:
        for e in p:
            if len(flip_edge(e,p,ppts)) != 1:
                unflippable[e] = unflippable.get(e,[]) + [p]
    return unflippable
