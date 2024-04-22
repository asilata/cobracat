# load("../complexes.sage")
# load("../zigzagalgebra.sage")
# load("../zigzagmodules.sage")
# load("../braidactions.sage")

import itertools
import networkx as nx

load("setup_cartan.sage")

ct = CartanType("D4")
p, s, burau = zz_setup(ct)

W = WeylGroup(ct)

reduced_words = sorted([x.reduced_word() for x in W], key=len)

def all_signed_combinations(w):
    """
    Return all possible signed lifts of a given word w in the Weyl group.
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

braid_lifts = sum([all_signed_combinations(w) for w in reduced_words], [])

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
    

def generate_elements_in_heart():
    admissible_braid_lifts = [tuple(x) for x in braid_lifts]
    outputs = []
    long_word_length = max([len(x) for x in reduced_words])
    for l in range(0,long_word_length + 1):
        lifts_of_length_l = [x for x in admissible_braid_lifts if len(x) == l]
        for b in lifts_of_length_l:
            b_of_p1 = composeAll([s[i] for i in b])(p[1])
            if not in_heart(b_of_p1):
                print(len(admissible_braid_lifts))
                admissible_braid_lifts = [x for x in admissible_braid_lifts if len(x) > l and x[-len(b):] != b]
                print(len(admissible_braid_lifts))                
            else:
               outputs = outputs + [b_of_p1]
    return outputs

def equal_upto_shift(x,y):
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

def de_duplicate(complexes_list):
    new_list = []
    for a in complexes_list:
        equal_elements = [n for n in new_list if equal_upto_shift(a,n)]
        if len(equal_elements) > 0:
            continue
        else:
            new_list = new_list + [a]
    return new_list

def exists_hom_one(A,B):
    h = hom(A,B).qPolynomial()
    hom_ones = [c for c in h.coefficients() if c[1] == -1]
    if len(hom_ones) >= 1:
        return True
    else:
        return False

def find_no_hom_one_cliques(candidates):
    no_hom_one_pairs = [(candidates[i], candidates[j]) for i in range(0, len(candidates)) for j in range(0,i) if not exists_hom_one(candidates[i], candidates[j])]
    G = nx.Graph()
    G.add_edges_from(no_hom_one_pairs)
    g_cliques = list(nx.find_cliques(G))
    m = max([len(x) for x in g_cliques])
    return g_cliques, [x for x in g_cliques if len(x) == m]

def find_external_edges(candidates):
    externals = []
    for x in candidates:
        hom_ones = [y for y in candidates if exists_hom_one(x,y)]
        if len(hom_ones) > 0:
            continue
        else:
            externals = externals + [x]
    return externals

def find_heart():
    heart = generate_elements_in_heart()
    smaller_heart_elements = de_duplicate(heart)
    return smaller_heart_elements
    
def find_ppts():
    smaller_heart_elements = find_heart()
    ppts = [frozenset(x) for x in find_no_hom_one_cliques(smaller_heart_elements)]
    return smaller_heart_elements, ppts

def flip_graph(ppts):
    max_intersection_pairs = [(ppts[i], ppts[j]) for i in range(0, len(ppts)) for j in range(0,i) if len(ppts[i] & ppts[j]) == len(ppts[i]) - 1]
    G = nx.Graph()
    G.add_edges_from(max_intersection_pairs)
    return G

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
