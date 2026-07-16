r"""
Semi-stable objects under given stability conditions

Given a root in the root system associated to some Cartan type, and a stability condition (enough to mention the central charge on the simple roots), we calculate the semi-stable object(s) under the given stability condition.

AUTHORS:

- Tanisha Talekar (2025-02-17): initial version

- Asilata Bapat (2027-06-16)

"""

# ****************************************************************************
#       Copyright (C) 2025 Tanisha Talekar <tanishatalekar23@iisertvm.ac.in> 
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

load("zigzagalgebra_element.pyx")
load("zigzagalgebra.sage")
load("projective-zigzagmodules.sage")
load("complexes.sage")
load("braidactions.sage")

import itertools

def root_sequence(root, with_writing=False):
    r"""
    Given a positive root, it calculates a root sequence corresponding to some writing of the root as the image under of a simple root under a reduced word.
    
    INPUT:
    
    `root` -- An element of a RootSpace
    `with_writing` -- a Boolean to indicate whether the corresponding writing should also be output. Default False.
    
    
    OUTPUT:
    
    A root sequence as a list of roots. If `with_writing` is true, then also output a tuple (index, word) such that the given root is obtained by applying `word` to the simple root at `index`.
    
    EXAMPLES:
    
        sage: alpha = RootSystem(['D',5]).root_space().simple_roots()
        sage: root_sequence(alpha[1])
        [alpha[1]]

        sage: root_sequence(alpha[1] + 2*alpha[2] + 2*alpha[3] + alpha[4] + alpha[5], with_writing=True)
        ([alpha[2],
        alpha[1] + alpha[2],
        alpha[2] + alpha[3],
        alpha[1] + alpha[2] + alpha[3],
        alpha[2] + alpha[3] + alpha[4],
        alpha[1] + alpha[2] + alpha[3] + alpha[4],
        alpha[1] + 2*alpha[2] + 2*alpha[3] + alpha[4] + alpha[5]],
        (5, (2, 1, 3, 2, 4, 3)))
    """
    root_space = root.parent()
    alpha = root_space.simple_roots()

    # Make sure the given root is a positive root.
    if root not in root_space.positive_roots():
        raise ValueError("{} is not a positive root of {}!".format(root, root_space))

    # Find a reduced_word that sends the simple root at index to the given root.
    index, reduced_word = root.to_simple_root(reduced_word=True)
    root_sequence = []
    for i in range(len(reduced_word)):
        root_sequence.append(alpha[reduced_word[i]].weyl_action(reduced_word[:i]))

    root_sequence.append(alpha[index].weyl_action(reduced_word))

    if with_writing:
        return root_sequence, (index, reduced_word)
    else:
        return root_sequence

class CentralCharge(object):
    """
    A class to encode the central charge associated to a given Cartan Type.
    """
    def __init__(self, ct, simple_values):
        self.cartan_type = CartanType(ct)
        self.root_space = RootSystem(ct).root_space()

        if len(simple_values) != len(self.cartan_type.index_set()):
            raise ValueError("The given sequence of simple values {} does not match the index set {} of the Cartan type {}!".format(simple_values, self.cartan_type.index_set(), self.cartan_type))
        self.simple_values = simple_values
        self.simple_roots = self.root_space.simple_roots()
        self.positive_roots = self.root_space.positive_roots()

    @cached_method
    def __getitem__(self, rs_element):
        if rs_element not in self.root_space:
            raise ValueError("{} is not an element of the domain {} of this central charge!".format(rs_element, self.root_space))
        
        rs_vector = rs_element.to_vector()
        return sum([a * b for (a,b) in zip(rs_vector, self.simple_values)])

    def __str__(self):
        return "Central charge on {} taking values {} on the simple roots".format(RootSystem(self.cartan_type), self.simple_values)

    def __repr__(self):
        return str(self)

    @cached_method
    def phase_order(self):
        """
        Return a list of tuples (positive_root, self[positive_root]) in order of increasing phase.
        """
        cc_positive_values = [(r, self[r]) for r in self.positive_roots]

        return sorted(cc_positive_values, key = lambda x : float(arg(x[1])))

def braid_lifts(root, cc):
    rs, (index, word) = root_sequence(root, with_writing=True)
    comparison_vector = cc[rs[-1]]
    signs = []
    for r in rs[:-1]:
        if arg(cc[r]) < arg(comparison_vector):
            signs.append([-1])
        elif arg(cc[r]) > arg(comparison_vector):
            signs.append([1])
        else:
            signs.append([1,-1])
    all_sign_vectors = itertools.product(*signs)
    return index, [list(zip(word,x)) for x in all_sign_vectors]
        
class StandardStabilityCondition(object):
    """
    A class to encode a standard stability condition of some CartanType.
    """
    def __init__(self, central_charge):
        self.central_charge = central_charge
        self.cartan_type = central_charge.cartan_type
        self.root_space = central_charge.root_space
        self.zzalgebra = ZigZagAlgebra(self.cartan_type)
        self.gen_projectives = {}
        for i in self.cartan_type.index_set():
            fpi = ProjectiveZigZagModule(self.zzalgebra, i)
            pi = ProjectiveComplex(self.zzalgebra)
            pi.add_object_at(0, fpi)
            self.gen_projectives[i] = pi

    @cached_method
    def __getitem__(self, root):
        """
        Return either the unique stable object of (positive) class root, or a list of semi-stable objects of class root.
        """
        if root not in self.root_space.positive_roots():
            raise ValueError("{} is not a positive root in {}!".format(root, self.root_space))
        if root in self.root_space.simple_roots():
            return [self.gen_projectives[root.to_simple_root()]]

        index, braid_words = braid_lifts(root, self.central_charge)
        return [self.gen_projectives[index].braid_action(w) for w in braid_words]

    def __str__(self):
        return "A standard stability condition on {} with {}".format(self.root_space.cartan_type(), self.central_charge)

    def __repr__(self):
        return str(self)

    @cached_method
    def phase_order(self):
        """
        Return a list of tuples (positive root, [semistables of that root class]) in order of increasing phase.
        """
        return [(p[0], self[p[0]]) for p in self.central_charge.phase_order()]

def find_top(obj, ssc):
    """
    Compute the top HN piece of the given object with respect to the stability condition being considered.

    INPUT:

    - `obj` : a ProjectiveComplex
    - `ssc`: a StandardStabilityCondition

    OUTPUT:

    A tuple (x, h) where x is a stable object of the highest phase that has a non-trivial morphism to obj.
    The morphism h is one such morphism.
    """
    stab = flatten([x[1] for x in ssc.phase_order()])
    
    highest_t_deg_found = -Infinity
    highest_object_found = None
    highest_hom_found = None

    for stable in stab:
        (h, hom_data)= stable.hom(obj, with_homs=True)
        h.minimize()

        for i in range(h.min_index, h.max_index+1):
            for j in range(0, len(h[i])):
                if h[i][j].graded_degree - i >= highest_t_deg_found:
                    highest_t_deg_found = h[i][j].graded_degree - i
                    highest_object_found = stable.graded_shift(h[i][j].graded_degree).homological_shift(-i)
                    highest_hom_found = {x+i : hom_data[i][j][x] for x in hom_data[i][j]}
                    
    return (highest_object_found, highest_hom_found)

def hn_factors(obj, ssc):
    """
    Compute the HN filtration of object obj with respect to the given (generic) stability condition.

    INPUT:

    - `obj` : a ProjectiveComplex.
    - `ssc`: a sequence of objects ordered by increasing phase

    OUTPUT:

    A list of stable HN (or refined HN) factors of obj, in non-increasing phase order.
    """
    hn_factors_list = []
    
    while not obj.is_zero():
        (top_piece, top_hom) = find_top(obj, ssc)
        hn_factors_list.append(top_piece)
        obj = cone(top_piece, obj, top_hom).homological_shift(1)
        obj.minimize()

    return hn_factors_list

