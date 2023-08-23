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
        Arguments -
        algebra: Algebra (not necessarily commutative) over which the elements of the complex are modules 
        objects: A dictionary of type {i: [P]} where i is an integer and P is projective module over the base ring. A projective module can be anything that implements the following methods: 
        (1) P.is_annihilated_by(r) : Is right multiplication by r the zero map on P?
        (2) P.is_invertible(r): Is right multiplication by r an invertible map on P?
        (3) P.invert(r): If right multiplication by r is invertible, return a ring element which acts as its inverse.
        (4) P.hom(Q): Returns a set of ring elements (monomials), which by right multiplication, define a basis of Hom(P,Q)
        See the class ProjectiveModuleOverField for an example.

        maps: A dictionary of type {i: {(a,b): r}} where i is an integer, (a,b) is a pair of positive integers, and r is an element of algebra. The key i represents the map from the i-th to (i+1)-th index of the complex. The pair (a,b) says that the map is from the a-th object in the list of objects at the i-th index to the b-th object in the list of objects at the (i+1)-th index. The value r says that the map is given by right multiplication by r. Currently, there is no provision to specify more complicated maps. 

        names: A dictionary of type {P: n} where P is a projective module and n is a string that acts as a short name for P to put in the string representation of the complex.
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
                s = s + "+".join([self._names[x] if x in self._names.keys() else str(x) for x in objects])
            if i < largest:
                s = s + " → "
        s = s + " :[" + str(largest) + "]"
        return s

    def __repr__(self):
        return str(self)

    def objects_at_index(self, i):
        r"""
        A list of projective objects whose direct sum forms the i-th component of self.
        """        
        return list(self.objects.get(i, []))

    def maps_at_index(self, i):
        r"""
        The map from the i-th to the (i+1)-th component, represented as a dictionary {(a,b): r}. The map corresponds to (right) multiplication by the matrix associated to this dictionary.
        """
        return self.maps.get(i, {}).copy()

    def names(self):
        r"""
        Short, readable names for the projective objects in the complex.
        """
        return self._names.copy()

    def homological_shift_by(self, n = 1):
        r"""
        A new complex obtained by homologically shifting self by [n].
        """
        new_objects = {x-n: self.objects[x] for x in self.objects.keys()}
        new_maps = {x-n: {k: (-1)^n * self.maps[x][k] for k in self.maps[x].keys()} for x in self.maps.keys()}
        return ProjectiveComplex(self.algebra, new_objects, new_maps, self._names)

    def copy(self):
        return ProjectiveComplex(self.algebra, self.objects, self.maps, self._names)

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
            self._names[obj] = name

        self.min_index = min(index, self.min_index)
        self.max_index = max(index, self.max_index)


    def q_polynomial(self, variables = lazy_list(var('q') for i in count())):
        answer = 0
        self.minimize()
        for i in range(self.min_index, self.max_index+1):
            rest_variables = lazy_list(variables[i+1] for i in count())
            answer = answer + variables[0]^(-i) * sum([obj.q_polynomial(rest_variables) for obj in self.objects_at_index(i)])
        return answer.expand()

    # def getLevels(self, i):
    #     '''
    #     Return a list of the possible levels of objects at index i.
    #     '''
    #     if i < self.min_index or i > self.max_index:
    #         return None
    #     return [ob.twist() - i for ob in self.objects_at_index(i)]

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
            obj = self.objects_at_index(i)
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
            for m in self.maps_at_index(i):
                D.add_edge((obj_names[(i,m[0])], obj_names[(i+1,m[1])]), label=str(self.maps_at_index(i)[m]))

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
        for i in list(self.maps.keys()):
            for k in list(self.maps[i].keys()):
                if self.maps[i][k] == 0:
                    self.maps[i].pop(k)
        # Remove objects
        for i in list(self.objects.keys()):
            if self.objects[i] == []:
                self.objects.pop(i)

        if len(self.objects.keys()) > 0:
            self.min_index = min(self.objects.keys())
            self.max_index = max(self.objects.keys())
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
        if i < 0 or i >= len(self.objects_at_index(index)):
            raise IndexError("Index out of bounds")
        if j < 0 or j >= len(self.objects_at_index(index+1)):
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
                if not self.objects_at_index(i)[k[0]].is_annihilated_by(v):
                    print("Differential squared not zero at " + str(i) + ".")
                return False

        return True

    def direct_sum(self, Q):
        '''
        Direct sum of this complex and the complex Q
        '''
        # By convention, the objects of Q go after the objects of self, in order.
        objs, maps = {}, {}
        names = self.names()
        names.update(Q.names())
        smallest = min([self.min_index,Q.min_index])
        largest = max([self.max_index,Q.max_index])

        for i in range(smallest, largest + 1):
            objs[i] = self.objects_at_index(i) + Q.objects_at_index(i)

        for k in range(smallest, largest):
            maps[k] = self.maps_at_index(k)
            l,w = len(self.objects_at_index(k)), len(self.objects_at_index(k+1))
            for (p,q) in Q.maps_at_index(k):
                maps[k][(p+l,q+w)] = Q.maps_at_index(k)[(p,q)]
        return ProjectiveComplex(self.algebra, objs, maps, names)
                      
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
            maps = self.maps_at_index(index)
            objects = self.objects_at_index(index)
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
            alphaInverse = self.objects_at_index(index)[source].invert(alpha)
            sourceBasis = self.objects_at_index(index)[source].basis()
            for i in range(0, len(self.objects_at_index(index))):
                for j in range(0, len(self.objects_at_index(index+1))):
                    # Update any maps from i to j by adding the correction factor, unless (i,j) = (source, target).
                    if (i,j) == (source, target):
                        changeij = 0
                    else:
                        changeij = self.maps_at_index(index).get((i,target), 0) * alphaInverse * self.maps_at_index(index).get((source,j), 0)
                    newMapsIndex[(i,j)] = self.maps_at_index(index).get((i,j), 0) - changeij
                    
                # For each object at index except for the source, update the basis if any
                # to reflect the splitting.
                if sourceBasis != None and i != source:
                    iBasis = self.objects_at_index(index)[i].basis()
                    if iBasis != None:
                        change = self.maps_at_index(index).get((i,target), 0) * alphaInverse
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
            matrixAtIndex = matrix(len(self.objects_at_index(index))+1, len(self.objects_at_index(index+1))+1, self.maps_at_index(index))
            newMatrixAtIndex = matrixAtIndex.delete_rows([source]).delete_columns([target])
            self.maps[index] = newMatrixAtIndex.dict()

            matrixAtIndexMinus1 = matrix(len(self.objects_at_index(index-1)), len(self.objects_at_index(index))+1, self.maps_at_index(index-1))
            if matrixAtIndexMinus1.ncols() > 0:
                newMatrixAtIndexMinus1 = matrixAtIndexMinus1.delete_columns([source])
                self.maps[index-1] = newMatrixAtIndexMinus1.dict()

            matrixAtIndexPlus1 = matrix(len(self.objects_at_index(index+1))+1, len(self.objects_at_index(index+2)) ,self.maps_at_index(index+1))
            if matrixAtIndexPlus1.nrows() > 0:
                newMatrixAtIndexPlus1 = matrixAtIndexPlus1.delete_rows([target])
                self.maps[index+1] = newMatrixAtIndexPlus1.dict()

        #Finally we do a cleanup
        self.cleanup()
        return 

# The hom complex
def hom(P, Q, degree=0):
    Z = P.algebra
    def inducedMap1(P1, P2, Q, r, sign=1):
        #The map from Hom(P2, Q) -> Hom(P1, Q) induced by r: P1 -> P2
        matrix = {}
        P2Q = P2.hom(Q)
        P1Q = P1.hom(Q)
        for i in range(0, len(P2Q)):
            h = P2Q[i]
            rh = (r*h).monomial_coefficients()
            for basis_element_index in rh.keys():
                if rh[basis_element_index] != 0:                                
                    basis_element = Z.basis()[basis_element_index]
                    if basis_element not in P1Q:
                        raise TypeError("Unrecognized hom: " + str(basis_element) +" not in " + P1Q)
                    j = P1Q.index(basis_element)
                    matrix[(i,j)] = rh[basis_element_index] * sign
        return matrix

    def inducedMap2(P, Q1, Q2, r, sign=1):
        #The map from Hom(P, Q1) -> Hom(P, Q2) induced by r: Q1 -> Q2
        matrix = {}
        PQ1 = P.hom(Q1)
        PQ2 = P.hom(Q2)
        for i in range(0, len(PQ1)):
            h = PQ1[i]
            hr = (h*r).monomial_coefficients()
            for basis_element_index in hr.keys():
                if hr[basis_element_index] != 0:                
                    basis_element = Z.basis()[basis_element_index]
                    if basis_element not in PQ2:
                        raise TypeError("Unrecognized hom: " + str(basis_element) +" not in " + str(PQ2))
                    j = PQ2.index(basis_element)
                    matrix[(i,j)] = hr[basis_element_index] * sign
        return matrix
    
    Q = Q.homological_shift_by(degree)

    doubleComplexObjects = {}
    for i in range(P.min_index, P.max_index+1):
        for j in range(Q.min_index, Q.max_index+1):
            doubleComplexObjects[(i,j)] = [(a,b) for a in range(0,len(P.objects_at_index(i))) for b in range(0,len(Q.objects_at_index(j)))]

    doubleComplexMaps = {}
    for (i,j) in doubleComplexObjects.keys():
        for (a,b) in doubleComplexObjects.get((i,j),[]):
            for (c,d) in doubleComplexObjects.get((i,j+1),[]):
                if a == c:
                    im2 = inducedMap2(P.objects_at_index(i)[a], Q.objects_at_index(j)[b], Q.objects_at_index(j+1)[d], Q.maps_at_index(j).get((b,d), 0))
                    doubleComplexMaps[((i,j,a,b),(i,j+1,c,d))] = im2

            for (c,d) in doubleComplexObjects.get((i-1,j),[]):
                if b == d:
                    im1 = inducedMap1(P.objects_at_index(i-1)[c], P.objects_at_index(i)[a], Q.objects_at_index(j)[b], P.maps_at_index(i-1).get((c,a),0), sign = -1 * (-1)**(i-j))
                    doubleComplexMaps[((i,j,a,b),(i-1,j,c,d))] = im1
                    
    # Collapsing the double complex to a single complex
    Z = P.algebra
    k = Z.base_ring()
    homComplex = ProjectiveComplex(algebra=k)
    renumberingDictionary = {}
    for (i,j) in doubleComplexObjects:
        for (a,b) in doubleComplexObjects[(i,j)]:
            Source = P.objects_at_index(i)[a]
            Target = Q.objects_at_index(j)[b]
            homs = Source.hom(Target)
            for hom_index in range(0, len(homs)):
                # Add a copy of the base field for each hom in the correct degree.
                # Store the hom as a basis, along with source and target indices.
                homComplex.add_object_at(j-i, GradedProjectiveModuleOverField(k,
                                                                          1,
                                                                          - Z.deg(homs[hom_index]) + Target.twist() - Source.twist(),
                                                                          name="k",
                                                                          basis={((i,a), (j,b)): homs[hom_index]}))
                # Remember where the object is stored.
                renumberingDictionary[(i,j,a,b,hom_index)] = len(homComplex.objects_at_index(j-i))-1 
    # Now add maps
    for ((i,j,a,b), (I,J,A,B)) in doubleComplexMaps.keys():
        maps = doubleComplexMaps[((i,j,a,b), (I,J,A,B))]
        for (alpha, beta) in maps.keys():
            homComplex.add_map_at(j-i, renumberingDictionary[(i,j,a,b,alpha)], renumberingDictionary[(I,J,A,B,beta)], maps[(alpha,beta)])
    return homComplex
            

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
            D.add_map_at(index, i, j+len(P.objects_at_index(index+1)), M[index][(i,j)])

    return D
    

def is_chain_map(P, Q, M):
    '''
    Check that M defines a map of chain complexes from P to Q.
    M must have the type {i: d} where i is an integer and d is a dictionary {(a,b): r} whose associated matrix defines the map from P to Q (by right multiplication).
    '''

    min_index = min(P.min_index, Q.min_index)
    max_index = max(P.max_index, Q.max_index)
    for i in range(min_index, max_index):
        dPi = matrix(len(P.objects_at_index(i)), len(P.objects_at_index(i+1)), P.maps_at_index(i))
        dQi = matrix(len(Q.objects_at_index(i)), len(Q.objects_at_index(i+1)), Q.maps_at_index(i))
        Mi = matrix(len(P.objects_at_index(i)), len(Q.objects_at_index(i)), M.get(i,{}))
        Mip1 = matrix(len(P.objects_at_index(i+1)), len(Q.objects_at_index(i+1)), M.get(i+1,{}))
        for k, v in (dPi*Mip1 - Mi*dQi).dict().items():
            if not P.objects_at_index(i)[k[0]].is_annihilated_by(v):
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

class GradedProjectiveModuleOverField(object):
    ''''
    The class of projective modules over a field (also known as graded vector spaces).
    '''
    # What is a basis?
    # Basis is supposed to be a dictionary {x:v} where the keys x could be anything, but the values v should be elements of a k-module.
    # The dictionary {x:v} is supposed to represent the element formal sum v*x.
    def __init__(self, basefield, dimension, grade = 0, name=None, basis=None):
        if not basefield.is_field():
            raise TypeError("Basefield not a field.")
        if not dimension.is_integral() or dimension < 0:
            raise TypeError("Invalid dimension.")
        self._vsp = VectorSpace(basefield,dimension)
        self._grade = grade
        if dimension > 1 and basis != None:
            raise NotImplementedError("Basis only implemented for one dimensional spaces")
        self._basis = basis
        if name:
            self._name = name
        else:
            self._name = self._vsp.__str__()

    def __str__(self):
        return self._name + "<" + str(self._grade) +">"

    def __repr__(self):
        return self._name + "<" + str(self._grade) +">"

    @cached_method
    def is_annihilated_by(self, r):
        return (r == 0)

    @cached_method
    def grade(self):
        return self._grade

    @cached_method
    def basis(self):
        return self._basis

    @cached_method    
    def is_invertible(self, r):
        return (r != 0)

    @cached_method
    def hom(self, Q):
        return [1]

    @cached_method    
    def q_polynomial(self, variables = [var('q')]):
        return variables[0]^self._grade

    @cached_method    
    def invert(self, r):
        if r != 0:
            return 1/r
        else:
            raise TypeError("Not invertible.")

