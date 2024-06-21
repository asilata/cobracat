r"""
The homotopy category of complexes of projective modules over an algebra

Given an algebra over a (usually commutative) base ring, we implement complexes of projective modules over that algebra. We assume that a projective module is anything that implements the following methods: 
        (1) P.is_annihilated_by(r) : Is right multiplication by r the zero map on P?
        (2) P.is_invertible(r): Is right multiplication by r an invertible map on P?
        (3) P.invert(r): If right multiplication by r is invertible, return a ring element which acts as its inverse.
        (4) P.hom(Q): Returns a set of ring elements (monomials), which by right multiplication, define a basis of Hom(P,Q)
See the classes `ProjectiveModuleOverField` and `GradedProjectiveModuleOverField` as examples.

AUTHORS:

- Asilata Bapat (2023-08-23): initial version

- Anand Deopurkar (2023-08-23): initial version

"""

# ****************************************************************************
#       Copyright (C) 2023 Asilata Bapat <asilata@alum.mit.edu> 
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

###################################################################
# The homotopy category of projective complexes over a fixed ring #
###################################################################
from sage.misc.lazy_list import lazy_list
from itertools import count

# Setting this to true calls is_chain_complex() each time a complex is created, as well as is_chain_map() when a cone is created.
DEBUG = False

class ProjectiveComplex(object):
    r"""
    The class of ProjectiveComplexes over a ring (not necessarily
    commutative). A projective complex is a complex of projective objects. 
    """
    def __init__(self, algebra, objects = {}, maps = {}, names={}):
        r"""
        INPUT:
        
        - `algebra` -- Algebra (not necessarily commutative) over which the elements of the complex are modules 
        - `objects` -- A dictionary of type `{i: [P]}` where `i` is an integer and `P` is projective module over the base ring.
            A projective module can be anything that implements the following methods: 
           1. `P.is_annihilated_by(r)` : Is right multiplication by `r` the zero map on `P`?
           2. `P.is_invertible(r)`: Is right multiplication by `r` an invertible map on `P`?
           3. `P.invert(r)`: If right multiplication by `r` is invertible, return a ring element which acts as its inverse.
           4. `P.hom(Q)`: Returns a set of ring elements (monomials), which by right multiplication, define a basis of `Hom(P,Q)`
        See the class `ProjectiveModuleOverField` for an example.

        - `maps` -- A dictionary of type `{i: {(a,b): r}}` where `i` is an integer, `(a,b)` is a pair of positive integers, and `r` is an element of algebra.
                    The key `i` represents the map from the `i`-th to `(i+1)`-th index of the complex.
                    The pair `(a,b)` says that the map is from the `a`-th object in the list of objects at the `i`-th index to the `b`-th object in the list of objects at the (i+1)-th index.
                    The value `r` says that the map is given by right multiplication by `r`.
                    Currently, there is no provision to specify more complicated maps. 

        - `names` -- A dictionary of type `{P: n}` where `P` is a projective module and `n` is a string that acts as a short name for `P` to put in the string representation of the complex.
        """
        self.objects = objects.copy()
        self.maps = maps.copy()
        self.names = names.copy()
        self.algebra = algebra
        if len(objects.keys()) > 0:
            self.min_index = min(objects.keys())
            self.max_index = max(objects.keys())
        else:
            self.min_index = 0
            self.max_index = 0

        if DEBUG:
            if not self.is_chain_complex():
                print("Warning: This is not a chain complex!")

    def __str__(self):
        ks = self.objects.keys()
        if len(ks) == 0:
            smallest,largest = 0,0
        else:
            smallest,largest = min(ks),max(ks)
        s = "[" + str(smallest) + "]: "

        for i in range(smallest,largest + 1):
            objects = self.objects.get(i,[])
            if len(objects) == 0:
                s = s + "0"
            else:
                s = s + "+".join([self.names[x] if x in self.names.keys() else str(x) for x in objects])
            if i < largest:
                s = s + " → "
        s = s + " :[" + str(largest) + "]"
        return s

    def __repr__(self):
        return str(self)

    def homological_shift_by(self, n = 1):
        r"""
        A new complex obtained by homologically shifting self by [n].
        """
        new_objects = {x-n: self.objects[x] for x in self.objects.keys()}
        # Sage thinks that (-1)^(-k) is a rational number :( so we decided to bit hack instead.
        # The expression (1 - (n % 2) * 2) evaluates to (-1)^n.
        new_maps = {x-n: {k: (1 - (n % 2) * 2) * self.maps[x][k] for k in self.maps[x].keys()} for x in self.maps.keys()}        
        return ProjectiveComplex(self.algebra, new_objects, new_maps, self.names)

    def copy(self):
        return ProjectiveComplex(self.algebra, self.objects, self.maps, self.names)

    def add_object_at(self, index, obj, name = None):
        r"""
        Add object `obj` at index `index` with name `name`. It goes as the last entry of the list of objects at index `index`.
        """
        if index not in self.objects:
            self.objects[index] = []
        self.objects[index].append(obj)

        if index not in self.maps:
            self.maps[index] = {}

        if name != None:
            self.names[obj] = name

        self.min_index = min(index, self.min_index)
        self.max_index = max(index, self.max_index)


    def q_polynomial(self, variables = lazy_list(var('q') for i in count())):
        answer = 0
        self.minimize()
        for i in range(self.min_index, self.max_index+1):
            rest_variables = lazy_list(variables[i+1] for i in count())
            answer = answer + variables[0]^(-i) * sum([obj.q_polynomial(rest_variables) for obj in self.objects.get(i,[])])
        return answer.expand()

    # def getLevels(self, i):
    #     '''
    #     Return a list of the possible levels of objects at index i.
    #     '''
    #     if i < self.min_index or i > self.max_index:
    #         return None
    #     return [ob.twist() - i for ob in self.objects.get(i,[])]

    # def minLevel(self):
    #     mins = [min(self.getLevels(i)) for i in range(self.min_index, self.max_index + 1)]
    #     return min(mins)

    # def maxLevel(self):
    #     maxs = [max(self.getLevels(i)) for i in range(self.min_index, self.max_index + 1)]
    #     return max(maxs)

    def show(self, **args):
        D = DiGraph()
        # Add objects
        obj_names = {}
        heights = {}
        level_sets = {}
        for i in range(self.min_index, self.max_index+1):
            heights[i] = []
            obj = self.objects.get(i,[])
            for j in range(0,len(obj)):
                ob = obj[j]
                obj_names[(i,j)] = "{0}({1},{2})".format(str(ob),i,j)
                D.add_vertex(obj_names[(i,j)])
                level = i - ob.twist()
                if level not in level_sets.keys():
                    level_sets[level] = []
                level_sets[level].append(obj_names[(i,j)])
                heights[i].append(obj_names[(i,j)])

        # Add edges
        for i in range(self.min_index, self.max_index):
            for m in self.maps.get(i,{}):
                D.add_edge((obj_names[(i,m[0])], obj_names[(i+1,m[1])]), label=str(self.maps.get(i,{})[m]))

        # Try to modify the positions to minimize intersections.
        plot = D.graphplot(heights=heights, layout="acyclic", save_pos=True)
        positions = D.get_pos()


        new_positions = {}
        # The algorithm to get new positions is a bit hacky, but the idea is the following.
        # First, we create add the reverses of all edges, effectively making the graph undirected.
        # Then we find traverse the longest (simple) path and add the vertices we encounter to a sequence.
        # At every height, we sort the vertices so that the ones encountered earlier are to the left.

        G = D.to_undirected().to_directed()
        walk_index = {}
        i = 0
        paths = G.all_simple_paths()
        paths.sort(key=len, reverse=true)
        for p in paths:
            for v in p:
                if v not in walk_index.keys():
                    walk_index[v] = i
                    i = i + 1
        
        # Traverse vertices that were not a part of a path.
        for v in G.vertices():
            if v not in walk_index.keys():
                walk_index[v] = i
                i = i + 1

        for vertices_at_h in heights.values():
            vertices_at_h.sort(key=lambda v: walk_index[v])
            positions_at_h = [positions[v] for v in vertices_at_h]
            positions_at_h.sort(key=lambda x: x[0])
            for i in range(0, len(vertices_at_h)):
                new_positions[vertices_at_h[i]] = positions_at_h[i]

        # Let us reverse (x,y) so that the complex shows up horizontally.
        flip = lambda x: (x[1], x[0])
        for k in new_positions.keys():
            new_positions[k] = flip(new_positions[k])

        D.set_pos(new_positions)
        # only choose the lighter colors
        colorlist = [c.rgb() for c in colors.values() if sum(c.rgb()) > 2.0] 
        vertex_colors = {colorlist[level]:level_sets[level] for level in level_sets.keys()}
        return D.plot(vertex_colors=vertex_colors, **args)
        
    def cleanup(self):
        '''
        Remove spurious matrix entries (zeros) and spurious object lists (empty lists).
        '''
        # Remove maps
        for i in self.maps:
            for k in list(self.maps[i].keys()):
                if self.maps[i][k] == 0:
                    self.maps[i].pop(k)
        # Remove objects
        for i in list(self.objects.keys()):
            if self.objects[i] == []:
                self.objects.pop(i)

        if len(self.objects) > 0:
            self.min_index = min(self.objects)
            self.max_index = max(self.objects)
        else:
            self.min_index = 0
            self.max_index = 0

    def minimize(self):
        '''
        Apply minimize_at(i) for all i.
        '''
        for i in range(self.min_index, self.max_index):
            self.minimize_at(i)

    def add_map_at(self, index, i, j, scalar):
        '''
        Add a map from the i-th object at index to the j-th object at index+1 given by right multiplication by scalar.
        '''
        # All actions are right actions!
        if i < 0 or i >= len(self.objects.get(index,[])):
            raise IndexError("Index out of bounds")
        if j < 0 or j >= len(self.objects.get(index+1,[])):
            raise IndexError("Index out of bounds")
        
        if index not in self.maps:
            self.maps[index] = {}
        self.maps[index][(i,j)] = self.algebra(scalar)

    def is_chain_complex(self):
        '''
        Check that this forms a chain complex.
        '''
        matrices = {}
        for i in range(self.min_index, self.max_index):
            source_dim = len(self.objects.get(i, []))
            target_dim = len(self.objects.get(i+1, []))
            matrices[i] = matrix(source_dim, target_dim, self.maps.get(i, {}))

        for i in range(self.min_index, self.max_index-1):
            for k, v in (matrices[i] * matrices[i+1]).dict().items():
                if not self.objects.get(i,[])[k[0]].is_annihilated_by(v):
                    print("Differential squared not zero at " + str(i) + ".")
                return False

        return True

    def direct_sum(self, Q):
        '''
        Direct sum of this complex and the complex Q
        '''
        # By convention, the objects of Q go after the objects of self, in order.
        objs, maps = {}, {}
        names = self.names
        names.update(Q.names)
        smallest = min([self.min_index,Q.min_index])
        largest = max([self.max_index,Q.max_index])

        for i in range(smallest, largest + 1):
            objs[i] = self.objects.get(i,[]) + Q.objects.get(i,[])

        for k in range(smallest, largest):
            maps[k] = self.maps.get(k,{}).copy()
            l,w = len(self.objects.get(k,[])), len(self.objects.get(k+1,[]))
            for (p,q) in Q.maps.get(k,{}):
                maps[k][(p+l,q+w)] = Q.maps.get(k,{})[(p,q)]
        return ProjectiveComplex(self.algebra, objs, maps, names)

    # TODO
    # There is an obscure bug in the following function. To reproduce:
    # Reset the copy() function of ProjectiveComplex to make a shallow copy.
    # Consider the object X = b(p[2]), where b = composeAll([s[1],s[2],s[2],s[2],s[1]]).
    # Set Y = X.copy().
    # Then do Y.minimize(). Inspect X to observe several zeros.
    # Try X.minimize_using_matrix()
    
    def minimize_using_matrix(self, with_qis=False):
        # We first transform the given complex into a giant matrix, whose rows and columns are indexed by all the objects in the given matrix.
        # We first store this correspondence.
        N = sum([len(x) for x in self.objects.values()])
        objects = [(k, i) for k in sorted(self.objects) for i in range(0,len(self.objects[k]))]
        objects_dict = dict(zip(range(0,N), objects))
        objects_reverse_dict = dict(zip(objects, range(0,N)))
        maps_dict = {(objects_reverse_dict[(i,a)], objects_reverse_dict[(i+1,b)]) : self.maps[i][(a,b)] for i in self.maps for (a,b) in self.maps[i]}
        M = matrix(N,N,maps_dict, sparse=True)

        # P will become the change of basis matrix later.
        # Why can't we define the identity matrix over the zigzag algebra?
        # Bug with error "ValueError: inconsistent number of rows: should be x but got y".
        # P = identity_matrix(self.algebra._base_ring, N, sparse=True) if with_qis else None
        if with_qis:
            P = matrix({(i,i):self.algebra(1) for i in range(N)}, sparse=True)
        else:
            P = None

        dropped_objects = set([])
        def _reduce_at(i,j,inverse):
            '''
            INPUT:

            - `i` -- a row coordinate of `M`
            - `j` -- a column coordinate of `M`
            - `inverse` -- the inverse (map) of M[i,j]

            OUTPUT:

            None.  Clear the i-th row and j-th column of M by Gaussian elimination.
            '''
            # For each l != i,j, replace column l (C_l) by C_l - C_j * inverse * M[i,l].
            # Explicitly, replace M[p,l] by M[p,l] - M[p,j] * inverse * M[i,l]
            # Also clear the jth row. (This is the the corresponding inverse row operation.)
            
            # For each p != i,j, replace row p (R_p) by R_p - M[p,j] * inverse * R_i
            # Explicitly, replace M[p,l] by M[p,l] - M[p,j] * inverse * M[i,l]
            # Also clear the ith column. (This is the the corresponding inverse column operation.)
            
            # Only do the calculations for those p for which M[p,j] is
            # non-zero and those l for which M[i,l] is non-zero.

            # Column operations.
            for l in M.nonzero_positions_in_row(i):
                weight_l = inverse * M[i,l]
                if l != j:                
                    for p in M.nonzero_positions_in_column(j):
                        M[p,l] = M[p,l] - M[p,j] * weight_l
                        
                    # Update P if required.
                    # Setting C[j,l] = - weight_l, we want to compute (I_N - C) * P.
                    # So replace P[j,-] by P[j,-] + \Sigma_l P[l,-] * weight_l
                    if P is not None:
                        for p in P.nonzero_positions_in_row(l):
                            P[j,p] = P[j,p] + P[l,p] * weight_l

            # Row operations.
            for p in M.nonzero_positions_in_column(j):
                weight_p = M[p,j] * inverse
                if p != i:                
                    for l in M.nonzero_positions_in_row(i):
                        M[p,l] = M[p,l] - weight_p * M[i,l]
                        
                    # Update P if required.
                    # Setting R[p,i] = - weight_p, we want to compute (I_N + R) * P.
                    # So, replace P[p,-] by P[p,-] - weight_p * P[i,-]
                    if P is not None:
                        for l in P.nonzero_positions_in_row(i):
                            P[p,l] = P[p,l] - weight_p * P[i,l]
            
            # Clear the ith column.
            for p in M.nonzero_positions_in_column(i):
                M[p,i] = 0

            # Clear the jth row.
            for l in M.nonzero_positions_in_row(j):
                M[j,l] = 0

            return

        for i in range(0,N):
            for j in M.nonzero_positions_in_row(i):
                source_object_index = objects_dict[i]
                source_object = self.objects.get(source_object_index[0],[])[source_object_index[1]]
                
                if source_object.is_invertible(M[i,j]):
                    inverse = source_object.invert(M[i,j])                    
                    _reduce_at(i,j,inverse) # TODO
                    dropped_objects.add(i)
                    dropped_objects.add(j)
                    break

        # We now create a new complex using the Gaussian eliminated M.
        # For that, we need to compute the new dictionary.
        new_objects_dict = {}

        # Assume that the entries of objects_dict are lexicographically sorted.
        homological_index = -Infinity
        for obj_index in range(0,N):
            obj = objects_dict[obj_index]

            if obj[0] > homological_index:
                drop_count = 0
                homological_index = obj[0]

            if obj_index in dropped_objects:
                drop_count += 1
                continue

            new_objects_dict[obj_index] = (obj[0], obj[1]-drop_count)

        reduced_complex = ProjectiveComplex(self.algebra, names=self.names)
        for obj_index in range(0,N):
            if obj_index not in dropped_objects:
                obj = objects_dict[obj_index]
                reduced_complex.add_object_at(obj[0], self.objects.get(obj[0],[])[obj[1]])
            
        for (i,j) in M.nonzero_positions():
            if not(i in dropped_objects or j in dropped_objects):
                new_source = new_objects_dict.get(i,None)
                new_target = new_objects_dict.get(j,None)
                assert (new_source and new_target)
                assert (new_target[0] == new_source[0] + 1)
                                
                reduced_complex.add_map_at(new_source[0], new_source[1], new_target[1], M[i,j])

        # Finally, clean up reduced_complex
        reduced_complex.cleanup()

        #Use generated P to create chain map which is a qis, and return it optionally.
        if P is not None:
            qis = {}
            for i in range(0,N):
                if i in dropped_objects:
                    continue
                new_obj = new_objects_dict.get(i,None)
                homological_index = new_obj[0]
                if homological_index not in qis:
                    qis[homological_index] = {}
                for j in P.nonzero_positions_in_row(i):
                    old_obj = objects_dict.get(j, None)
                    assert homological_index == old_obj[0]
                    qis[homological_index][(new_obj[1], old_obj[1])] = P[i,j]
            return reduced_complex, qis
        else:
            return reduced_complex
        
                      
    def minimize_at(self, index):
        '''
        Factor out a complex (presumably the biggest such) that is chain homotopic to zero and is concentrated in degrees index and index+1.
        '''
        k = self.algebra.base_ring()
        
        # Find an object at index and an object at (index+1) with an isomorphism between them.
        # The return value is (i, j, fij), where:
        # i is the index of the source object at index,
        # j is the index of the target object at (index + 1), and
        # fij is the isomorphism between them.
        def _findIso(index):
            maps = self.maps.get(index,{})
            objects = self.objects.get(index,[])
            zero = self.algebra(0)
            for (i,j) in maps:
                fij = maps.get((i,j), zero)
                if objects[i].is_invertible(fij):
                    return i, j, fij
            return None, None, None

        is_already_minimized = False
        while not is_already_minimized:
            source, target, alpha = _findIso(index)
            if source == None or target == None or alpha == None:
                #print("Nothing left to minimize at " + str(index))
                is_already_minimized = True
                continue

            # Change the maps from index to index+1
            newMapsIndex = {}
            alphaInverse = self.objects.get(index,[])[source].invert(alpha)
            sourceBasis = self.objects.get(index,[])[source].basis()
            for i in range(0, len(self.objects.get(index,[]))):
                for j in range(0, len(self.objects.get(index+1,[]))):
                    # Update any maps from i to j by adding the correction factor, unless (i,j) = (source, target).
                    if (i,j) == (source, target):
                        changeij = 0
                    else:
                        changeij = self.maps.get(index,{}).get((i,target), 0) * alphaInverse * self.maps.get(index,{}).get((source,j), 0)
                    newMapsIndex[(i,j)] = self.maps.get(index,{}).get((i,j), 0) - changeij
                    
                # For each object at index except for the source, update the basis if any
                # to reflect the splitting.
                if sourceBasis != None and i != source:
                    iBasis = self.objects.get(index,[])[i].basis()
                    if iBasis != None:
                        change = self.maps.get(index,{}).get((i,target), 0) * alphaInverse
                        for elt in sourceBasis:
                            iBasis[elt] = iBasis.get(elt,0) - change * sourceBasis[elt]

            # The maps from index-1 to index and index+1 to index+2 do not need to be changed substantially, apart from the indexing.
            # Now we update the maps
            self.maps[index] = newMapsIndex

            # At this point, our complex is a direct sum of F (source) -> F (target)
            # and another complex C.
            # We simply drop the source and the target
            self.objects[index].pop(source)
            self.objects[index+1].pop(target)


            # We now re-index as needed
            matrixAtIndex = matrix(len(self.objects.get(index,[]))+1, len(self.objects.get(index+1,[]))+1, self.maps.get(index,{}))
            newMatrixAtIndex = matrixAtIndex.delete_rows([source]).delete_columns([target])
            self.maps[index] = newMatrixAtIndex.dict()

            matrixAtIndexMinus1 = matrix(len(self.objects.get(index-1,[])), len(self.objects.get(index,[]))+1, self.maps.get(index-1,{}))
            if matrixAtIndexMinus1.ncols() > 0:
                newMatrixAtIndexMinus1 = matrixAtIndexMinus1.delete_columns([source])
                self.maps[index-1] = newMatrixAtIndexMinus1.dict()

            matrixAtIndexPlus1 = matrix(len(self.objects.get(index+1,[]))+1, len(self.objects.get(index+2,[])) ,self.maps.get(index+1,{}))
            if matrixAtIndexPlus1.nrows() > 0:
                newMatrixAtIndexPlus1 = matrixAtIndexPlus1.delete_rows([target])
                self.maps[index+1] = newMatrixAtIndexPlus1.dict()

        #Finally we do a cleanup
        self.cleanup()
        return 

    def hom(self, Q, degree = 0, name="k", with_homs=False):
        r"""
        The complex `hom(self, Q)` of degree `degree` (default: 0).

        Return the minimized hom complex `H` of maps from `P` to `Q`.
        The object at the ith place of `H` is a direct sum of shifted copies of k.
        Each copy of `k` represents a chain map from `P` to `Q`.
        
        If `with_homs` is true, then we also return `B`, a dictionary
        indexed by objects of `H`. The ith entry of B is a dictionary
        indexed by all objects at position i of `H`. The element
        B[i][a] is the chain map represented by shifted copy of k that
        lives at H.objects[i][a].
        
        INPUT:

        - `Q` -- an object of class `ProjectiveComplex`
        - `degree` -- an integer
        - `name` -- a name for the base field
        - `with_homs` -- a boolean, signifying whether we want to
          return the actual chain maps from `self` to Q.
        
        OUTPUT:

        Either `(H,B)` (if `with_homs` is True), or `H` (if `with_homs` is False).
        - `H` -- a `ProjectiveComplex` of graded modules over `k`
        - `B` -- a dictionary `{i:{b: M}}` where `M` is a chain map
          from `self` tensor `H[i][b]` to `Q` shifted by `degree`.
        """
        Z = self.algebra
        Z_basis = list(Z.basis())
        
        Q = Q.homological_shift_by(degree)

        # In all explanations that follow, let P = self.

        #######################################
        #    Objects of the double complex    # 
        #######################################
        
        # Its entry at (i,j) is Hom(P[i],Q[j]) = + Hom(P[i][a], Q[j][b])
        # Hom(P[i][a],Q[j][b]) is in turn a direct sum of (monomial) maps, which are indexed by elements of the base ring.
        # So double_complex_objects[(i,j)] is a list containing entries of the form (a, b, index_m), where a and b are as
        # above, and index_m is the index of a monomial map from P[i][a] to Q[j][b].

        double_complex_objects = {}
        for i in self.objects:
            pi = self.objects[i]
            for j in Q.objects:
                qj = Q.objects[j]
                double_complex_objects[(i,j)] = [(a,b,Z_basis.index(m)) for a in range(len(pi)) for b in range(len(qj)) for m in pi[a].hom(qj[b])]

        ##############################
        # Maps of the double complex #
        ##############################

        # There are two kinds of differentials: horizontal and vertical (using row,column syntax like a matrix).
        
        # For the horizontal one: (i,j) --> (i,j+1)
        # The object at (i,j) is Hom(P[i], Q[j]) = [(a,b,index_m)], a list ranging over all monomial maps m: P[i][a] --> Q[j][b].
        # Consider induced_map, which is the composition P[i][a] --> Q[j][b] --> Q[j+1][d] where the first map is m and the second map is the differential of Q.
        
        # The object at (i,j+1) is Hom(P[i], Q[j+1]) = [(c,d,index_m2)], a list ranging over all monomial maps P[i][c] --> Q[j+1][d].
        # For a map m2: P[i][a] --> Q[j+1][d], let coeff2 be the coefficient of m2 in induced_map.
        # Therefore set double_complex_maps_horizontal[(i,j)][((a,b,index_m), (a, d, index_m2))] = coeff_m2.
        double_complex_maps_horizontal = {}

        for (i,j) in double_complex_objects:
            # TODO unsure if this is the best place to add this entry in.
            if (i,j) not in double_complex_maps_horizontal:
                double_complex_maps_horizontal[(i,j)] = {}
                
            for (a,b,index_m) in double_complex_objects[(i,j)]:
                source = self.objects[i][a]
                for (c,d) in Q.maps.get(j,{}):
                    if c != b:
                        continue
                    differential = Q.maps[j][(b,d)]
                    target = Q.objects[j+1][d]
                    induced_map = source.idempotent *  Z_basis[index_m] * differential * target.idempotent # Goes from P[i][a] -> Q[j+1][d]
                    
                    # Write the induced map in terms of the basis of P[i][a].hom(Q[j+1][d])
                    induced_map_monomial_coeffs = induced_map.monomial_coefficients()

                    for m2 in source.hom(target):
                        index_m2 = Z_basis.index(m2)
                        coeffm2 = induced_map_monomial_coeffs[index_m2]
                        if coeffm2 != 0:
                            double_complex_maps_horizontal[(i,j)][((a,b,index_m), (a,d,index_m2))] = coeffm2

                
        # For the vertical one: (i,j) --> (i-1,j)
        # The object at (i,j) is Hom(P[i], Q[j]) = [(a,b,index_m)], a list ranging over all monomial maps m: P[i][a] --> Q[j][b].
        # Consider induced_map, which is the composition P[i-1][c] --> P[i][a] --> Q[j][b] where the first map is a sign times the differential of P, and the second map is m.
        # The sign is (-1)^(i-j-1).
        
        # The object at (i-1,j) is Hom(P[i-1], Q[j]) = [(c,d,index_m2)], a list ranging over all monomial maps P[i-1][c] --> Q[j][d]
        # For a map m2: P[i-1][c] --> Q[j][b], let coeff2 be the coefficient of m2 in induced_map.
        # Therefore set double_complex_maps_vertical[(i,j)][((a,b,index_m), (c,b,index_m2))] = coeff_m2

        double_complex_maps_vertical = {}

        for (i,j) in double_complex_objects:
            if (i,j) not in double_complex_maps_vertical:
                double_complex_maps_vertical[(i,j)] = {}
                
            for (a,b,index_m) in double_complex_objects[(i,j)]:
                target = Q.objects[j][b]
                for (c,d) in self.maps.get(i-1,{}):
                    if d != a:
                        continue
                    differential = self.maps[i-1][(c,a)]
                    source = self.objects[i-1][c]
                    induced_map = (-1)**((i-j) % 2) * source.idempotent * differential * Z_basis[index_m]  * target.idempotent # Goes from P[i-1][c] -> Q[j][b]
                    
                    # Write the induced map in terms of the basis of P[i-1,c].hom(Q[j][b])
                    induced_map_monomial_coeffs = induced_map.monomial_coefficients()

                    for m2 in source.hom(target):
                        index_m2 = Z_basis.index(m2)                        
                        coeffm2 = induced_map_monomial_coeffs[index_m2]
                        if coeffm2 != 0:
                            double_complex_maps_vertical[(i,j)][(a,b,index_m), (c,b,index_m2)] = coeffm2


        ###############################
        # Flatten to a single complex #
        ###############################

        # Add objects to the flattened complex.

        k = Z.base_ring()
        hom_complex = ProjectiveComplex(algebra=k)
        renumbering_dictionary = {}
        for (i,j) in double_complex_objects:
            if j-i not in renumbering_dictionary:
                renumbering_dictionary[j-i] = []
            
            for (a,b,index_m) in double_complex_objects[i,j]:
                renumbering_dictionary[j-i].append((i,j,a,b,index_m))
                # Add to hom_complex.objects[j-i] a graded
                # one-dimensional k-vector space k<d> corresponding to
                # the map m: P[i][a] -> Q[j][b], where
                # d = Q[j][b].graded_degree - P[i][a].graded_degree - Z.deg(m)
                d = Q.objects[j][b].graded_degree - self.objects[i][a].graded_degree - Z.deg(Z_basis[index_m])
                # TODO: Using Z.deg(...) is a hack.  Not all algebras have a .deg method.
                # Ideally, P.hom(Q) should give a graded basis, and d should come from there.
                k_d = GradedProjectiveModuleOverField(k, 1, d, name=name)
                hom_complex.add_object_at(j-i, k_d)

        # Add horizontal maps
        for (i,j) in double_complex_maps_horizontal:
            dcmh_ij = double_complex_maps_horizontal[(i,j)]
            for ((a,b,index_m),(c,d,index_m2)) in dcmh_ij:
                assert a == c
                # Insert a map at hom_complex.maps[j-i] from the right
                # index to the right index.
                index_source = renumbering_dictionary[j-i].index((i,j,a,b,index_m))
                index_target = renumbering_dictionary[j+1-i].index((i,j+1,c,d,index_m2))
                hom_complex.add_map_at(j-i, index_source, index_target, dcmh_ij[(a,b,index_m), (c,d,index_m2)])

        # Add vertical maps
        for (i,j) in double_complex_maps_vertical:
            dcmv_ij = double_complex_maps_vertical[(i,j)]
            for ((a,b,index_m),(c,d,index_m2)) in dcmv_ij:
                assert b == d
                # Insert a map at hom_complex.maps[j-i] from the right
                # index to the right index.
                index_source = renumbering_dictionary[j-i].index((i,j,a,b,index_m))
                index_target = renumbering_dictionary[j+1-i].index((i-1,j,c,d,index_m2))
                hom_complex.add_map_at(j-i, index_source, index_target, dcmv_ij[(a,b,index_m), (c,d,index_m2)])
                
        #
        if with_homs:
            B = {}
            hom_complex_reduced, qis = hom_complex.minimize_using_matrix(with_qis=True)
            for j_minus_i in hom_complex_reduced.objects:
                # B[i] is a dictionary of type
                # {index_reduced_hom: chain-map from self to Q represented by reduced_hom} 

                B[j_minus_i] = {}
                for (index_reduced_hom, index_original_hom) in qis[j_minus_i]:
                    if index_reduced_hom not in B[j_minus_i]:
                        B[j_minus_i][index_reduced_hom] = {}

                    reduced_hom = B[j_minus_i][index_reduced_hom]
                    coefficient = qis[j_minus_i][(index_reduced_hom, index_original_hom)]
                    (i,j,a,b,index_m) = renumbering_dictionary[j_minus_i][index_original_hom]
                    assert j_minus_i == j - i

                    if i not in reduced_hom:
                        reduced_hom[i] = {}

                    reduced_hom[i][(a,b)] = reduced_hom[i].get((a,b),0) + coefficient * Z_basis[index_m]
            return hom_complex_reduced, B
        else:
            return hom_complex.minimize_using_matrix(with_qis = False)
    
def cone(P, Q, M):
    '''
    The cone of M: P -> Q. 
    M must define a map of chain complexes from P to Q.
    '''
    if DEBUG:
        if not is_chain_map(P, Q, M):
            raise TypeError("Not a chain map. Cannot make a cone.")

    D = P.direct_sum(Q.homological_shift_by(-1))
    for index in M.keys():
        for (i,j) in M.get(index,{}):
            D.add_map_at(index, i, j+len(P.objects.get(index+1,[])), M[index][(i,j)])

    return D
    

def is_chain_map(P, Q, M):
    '''
    Check that M defines a map of chain complexes from P to Q.
    M must have the type {i: d} where i is an integer and d is a dictionary {(a,b): r} whose associated matrix defines the map from P to Q (by right multiplication).
    '''

    min_index = min(P.min_index, Q.min_index)
    max_index = max(P.max_index, Q.max_index)
    for i in range(min_index, max_index):
        dPi = matrix(len(P.objects.get(i,[])), len(P.objects.get(i+1,[])), P.maps.get(i,{}))
        dQi = matrix(len(Q.objects.get(i,[])), len(Q.objects.get(i+1,[])), Q.maps.get(i,{}))
        Mi = matrix(len(P.objects.get(i,[])), len(Q.objects.get(i,[])), M.get(i,{}))
        Mip1 = matrix(len(P.objects.get(i+1,[])), len(Q.objects.get(i+1,[])), M.get(i+1,{}))
        for k, v in (dPi*Mip1 - Mi*dQi).dict().items():
            if not P.objects.get(i,[])[k[0]].is_annihilated_by(v):
                return False
    return True
    
class ProjectiveModuleOverField(object):
    ''''
    The class of projective modules over a field (also known as vector spaces).
    '''
    def __init__(self, basefield, dimension):
        if not basefield.is_field():
            raise TypeError("Basefield not a field.")
        if not dimension.is_integral() or dimension < 0:
            raise TypeError("Invalid dimension.")
        self._vsp = VectorSpace(basefield,dimension)

    def __str__(self):
        return self._vsp.__str__()

    def __repr__(self):
        return self._vsp.__repr__()

    @cached_method
    def is_annihilated_by(self, r):
        return (r == 0)

    @cached_method
    def is_invertible(self, r):
        return (r != 0)

    @cached_method
    def hom(self, Q):
        return [1]

    @cached_method
    def invert(self, r):
        if r != 0:
            return 1/r
        else:
            raise TypeError("Not invertible.")
