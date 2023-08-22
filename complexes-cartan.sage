###################################################################
# The homotopy category of projective complexes over a fixed ring #
###################################################################
from sage.misc.lazy_list import lazy_list
from itertools import count

# Setting this to true calls is_chain_complex() each time a complex is created, as well as checkMap() when a cone is created.
DEBUG = False

class ProjectiveComplex(object):
    r"""
    The class of ProjectiveComplexes over a ring (not necessarily
    commutative). A projective complex is a complex of projective objects. 
    """
    def __init__(self, base_ring, objects = {}, maps = {}, names={}):
        r"""
        Arguments -
        base_ring: Base ring (not necessarily commutative)

        objects: A dictionary of type {i: [P]} where i is an integer and P is projective module over the base ring. A projective module can be anything that implements the following methods: 
        (1) P.is_zero(r) : Is right multiplication by r the zero map on P?
        (2) P.is_invertible(r): Is right multiplication by r an invertible map on P?
        (3) P.invert(r): If right multiplication by r is invertible, return a ring element which acts as its inverse.
        (4) P.hom(Q): Returns a set of ring elements (monomials), which by right multiplication, define a basis of Hom(P,Q)
        See the class ProjectiveModuleOverField for an example.

        maps: A dictionary of type {i: {(a,b): r}} where i is an integer, (a,b) is a pair of positive integers, and r is an element of base_ring. The key i represents the map from the i-th to (i+1)-th place of the complex. The pair (a,b) says that the map is from the a-th object in the list of objects at the i-th place to the b-th object in the list of objects at the (i+1)-th place. The value r says that the map is given by right multiplication by r. Currently, there is no provision to specify more complicated maps. 

        names: A dictionary of type {P: n} where P is a projective module and n is a string that acts as a short name for P to put in the string representation of the complex.
        """
        self.objects = objects.copy()
        self.maps = maps.copy()
        self.names = names.copy()
        self.base_ring = base_ring
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
        ks = self._objects.keys()
        if len(ks) == 0:
            smallest,largest = 0,0
        else:
            smallest,largest = min(ks),max(ks)
        s = "[" + str(smallest) + "]: "

        for i in range(smallest,largest + 1):
            objects = self._objects.get(i,[])
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

    def min_index(self):
        '''
        An integer n such that self.objects(m) = [] for m < n.
        '''
        return self._min_index

    def max_index(self):
        '''
        An integer n such that self.objects(m) = [] for m > n.
        '''
        return self._max_index

    def objects(self, i):
        '''
        A list of projective objects whose direct sum forms the i-th component of self.
        '''
        return list(self._objects.get(i, []))

    def maps(self, i):
        '''
        The map from the i-th to the (i+1)-th component, represented as a dictionary {(a,b): r}. The map corresponds to (right) multiplication by the matrix associated to this dictionary.
        '''
        return self._maps.get(i, {}).copy()

    def names(self):
        '''
        Short, readable names for the projective objects in the complex.
        '''
        return self._names.copy()

    def shift(self, n = 1):
        '''
        A new complex obtained by homologically shifting self by [n].
        '''
        return ProjectiveComplex(self._base_ring,
                                 {x-n: self._objects[x] for x in self._objects.keys()},
                                 {x-n: {k: (-1)^n * self._maps[x][k] for k in self._maps[x].keys()}
                                  for x in self._maps.keys()},
                                 self._names)

    def copy(self):
        return ProjectiveComplex(self._base_ring, self._objects, self._maps, self._names)

    def addObject(self, place, obj, name = None):
        '''
        Add object `obj` at place `place` with name `name`. It goes as the last entry of the list of objects at place `place`.
        '''
        if place not in self._objects:
            self._objects[place] = []
        self._objects[place].append(obj)

        if place not in self._maps:
            self._maps[place] = {}

        if name != None:
            self._names[obj] = name

        self._min_index = min(place, self._min_index)
        self._max_index = max(place, self._max_index)


    def qPolynomial(self, variables = lazy_list(var('q') for i in count())):
        answer = 0
        self.minimize()
        for i in range(self.min_index(), self.max_index()+1):
            restVariables = lazy_list(variables[i+1] for i in count())
            answer = answer + variables[0]^(-i) * sum([obj.qPolynomial(restVariables) for obj in self.objects(i)])
        return answer.expand()

    def getLevels(self, i):
        '''
        Return a list of the possible levels of objects at index i.
        '''
        if i < self.min_index() or i > self.max_index():
            return None
        return [ob.twist() - i for ob in self.objects(i)]

    def minLevel(self):
        mins = [min(self.getLevels(i)) for i in range(self.min_index(), self.max_index() + 1)]
        return min(mins)

    def maxLevel(self):
        maxs = [max(self.getLevels(i)) for i in range(self.min_index(), self.max_index() + 1)]
        return max(maxs)

    def show(self, **args):
        D = DiGraph()
        # Add objects
        objNames = {}
        heights = {}
        levelSets = {}
        for i in range(self.min_index(), self.max_index()+1):
            heights[i] = []
            obj = self.objects(i)
            for j in range(0,len(obj)):
                ob = obj[j]
                objNames[(i,j)] = "{0}({1},{2})".format(str(ob),i,j)
                D.add_vertex(objNames[(i,j)])
                level = i - ob.twist()
                if level not in levelSets.keys():
                    levelSets[level] = []
                levelSets[level].append(objNames[(i,j)])
                heights[i].append(objNames[(i,j)])

        # Add edges
        for i in range(self.min_index(), self.max_index()):
            for m in self.maps(i):
                D.add_edge((objNames[(i,m[0])], objNames[(i+1,m[1])]), label=str(self.maps(i)[m]))

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
        vertex_colors = {colorlist[level]:levelSets[level] for level in levelSets.keys()}
        return D.plot(vertex_colors=vertex_colors, **args)
        
    def cleanUp(self):
        '''
        Remove spurious matrix entries (zeros) and spurious object lists (empty lists).
        '''
        # Remove maps
        for i in list(self._maps.keys()):
            for k in list(self._maps[i].keys()):
                if self._maps[i][k] == 0:
                    self._maps[i].pop(k)
        # Remove objects
        for i in list(self._objects.keys()):
            if self._objects[i] == []:
                self._objects.pop(i)

        if len(self._objects.keys()) > 0:
            self._min_index = min(self._objects.keys())
            self._max_index = max(self._objects.keys())
        else:
            self._min_index = 0
            self._max_index = 0

    def minimize(self):
        '''
        Apply minimizeAt(i) for all i.
        '''
        for i in range(self.min_index(), self.max_index()):
            self.minimizeAt(i)
        

    def addMap(self, place, i, j, scalar):
        '''
        Add a map from the i-th object at place to the j-th object at place+1 given by right multiplication by scalar.
        '''
        # All actions are right actions!
        if i < 0 or i >= len(self.objects(place)):
            raise IndexError("Index out of bounds")
        if j < 0 or j >= len(self.objects(place+1)):
            raise IndexError("Index out of bounds")
        
        if place not in self._maps:
            self._maps[place] = {}
        self._maps[place][(i,j)] = self._base_ring(scalar)

    def is_chain_complex(self):
        '''
        Check that this forms a chain complex.
        '''
        matrices = {}
        for i in range(self.min_index(), self.max_index()):
            sourceDim = len(self._objects.get(i, []))
            targetDim = len(self._objects.get(i+1, []))
            matrices[i] = matrix(sourceDim, targetDim, self._maps.get(i, {}))

        for i in range(self.min_index(), self.max_index()-1):
            for k, v in (matrices[i] * matrices[i+1]).dict().items():
                if not self.objects(i)[k[0]].is_zero(v):
                    print("Differential squared not zero at " + str(i) + ".")
                return False

        return True

    def directSum(self, Q):
        '''
        Direct sum of this complex and the complex Q
        '''
        # By convention, the objects of Q go after the objects of self, in order.
        objs, maps = {}, {}
        names = self.names()
        names.update(Q.names())
        smallest = min([self.min_index(),Q.min_index()])
        largest = max([self.max_index(),Q.max_index()])

        for i in range(smallest, largest + 1):
            objs[i] = self.objects(i) + Q.objects(i)

        for k in range(smallest, largest):
            maps[k] = self.maps(k)
            l,w = len(self.objects(k)), len(self.objects(k+1))
            for (p,q) in Q.maps(k):
                maps[k][(p+l,q+w)] = Q.maps(k)[(p,q)]
        return ProjectiveComplex(self._base_ring, objs, maps, names)
                      
    def minimizeAt(self, place):
        '''
        Factor out a complex (presumably the biggest such) that is chain homotopic to zero and is concentrated in degrees place and place+1.
        '''
        k = self._base_ring.base_ring()
        
        # Find an object at place and an object at (place+1) with an isomorphism between them.
        # The return value is (i, j, fij), where:
        # i is the index of the source object at place,
        # j is the index of the target object at (place + 1), and
        # fij is the isomorphism between them.
        def _findIso(place):
            maps = self.maps(place)
            objects = self.objects(place)
            zero = self._base_ring(0)
            for (i,j) in maps:
                fij = maps.get((i,j), zero)
                if objects[i].is_invertible(fij):
                    return i, j, fij
            return None, None, None

        alreadyMinimized = False
        while not alreadyMinimized:
            source, target, alpha = _findIso(place)
            if source == None or target == None or alpha == None:
                #print("Nothing left to minimize at " + str(place))
                alreadyMinimized = True
                continue

            # Change the maps from place to place+1
            newMapsPlace = {}
            alphaInverse = self.objects(place)[source].invert(alpha)
            sourceBasis = self.objects(place)[source].basis()
            for i in range(0, len(self.objects(place))):
                for j in range(0, len(self.objects(place+1))):
                    # Update any maps from i to j by adding the correction factor, unless (i,j) = (source, target).
                    if (i,j) == (source, target):
                        changeij = 0
                    else:
                        changeij = self.maps(place).get((i,target), 0) * alphaInverse * self.maps(place).get((source,j), 0)
                    newMapsPlace[(i,j)] = self.maps(place).get((i,j), 0) - changeij
                    
                # For each object at place except for the source, update the basis if any
                # to reflect the splitting.
                if sourceBasis != None and i != source:
                    iBasis = self.objects(place)[i].basis()
                    if iBasis != None:
                        change = self.maps(place).get((i,target), 0) * alphaInverse
                        for elt in sourceBasis:
                            iBasis[elt] = iBasis.get(elt,0) - change * sourceBasis[elt]

            # The maps from place-1 to place and place+1 to place+2 do not need to be changed substantially, apart from the indexing.
            # Now we update the maps
            self._maps[place] = newMapsPlace

            # At this point, our complex is a direct sum of F (source) -> F (target)
            # and another complex C.
            # We simply drop the source and the target
            self._objects[place].pop(source)
            self._objects[place+1].pop(target)


            # We now re-index as needed
            matrixAtPlace = matrix(len(self.objects(place))+1, len(self.objects(place+1))+1, self.maps(place))
            newMatrixAtPlace = matrixAtPlace.delete_rows([source]).delete_columns([target])
            self._maps[place] = newMatrixAtPlace.dict()

            matrixAtPlaceMinus1 = matrix(len(self.objects(place-1)), len(self.objects(place))+1, self.maps(place-1))
            if matrixAtPlaceMinus1.ncols() > 0:
                newMatrixAtPlaceMinus1 = matrixAtPlaceMinus1.delete_columns([source])
                self._maps[place-1] = newMatrixAtPlaceMinus1.dict()

            matrixAtPlacePlus1 = matrix(len(self.objects(place+1))+1, len(self.objects(place+2)) ,self.maps(place+1))
            if matrixAtPlacePlus1.nrows() > 0:
                newMatrixAtPlacePlus1 = matrixAtPlacePlus1.delete_rows([target])
                self._maps[place+1] = newMatrixAtPlacePlus1.dict()

        #Finally we do a cleanup
        self.cleanUp()
        return 

# The hom complex
def hom(P, Q, degree=0):
    Z = P.base_ring()
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
    
    Q = Q.shift(degree)

    doubleComplexObjects = {}
    for i in range(P.min_index(), P.max_index()+1):
        for j in range(Q.min_index(), Q.max_index()+1):
            doubleComplexObjects[(i,j)] = [(a,b) for a in range(0,len(P.objects(i))) for b in range(0,len(Q.objects(j)))]

    doubleComplexMaps = {}
    for (i,j) in doubleComplexObjects.keys():
        for (a,b) in doubleComplexObjects.get((i,j),[]):
            for (c,d) in doubleComplexObjects.get((i,j+1),[]):
                if a == c:
                    im2 = inducedMap2(P.objects(i)[a], Q.objects(j)[b], Q.objects(j+1)[d], Q.maps(j).get((b,d), 0))
                    doubleComplexMaps[((i,j,a,b),(i,j+1,c,d))] = im2

            for (c,d) in doubleComplexObjects.get((i-1,j),[]):
                if b == d:
                    im1 = inducedMap1(P.objects(i-1)[c], P.objects(i)[a], Q.objects(j)[b], P.maps(i-1).get((c,a),0), sign = -1 * (-1)**(i-j))
                    doubleComplexMaps[((i,j,a,b),(i-1,j,c,d))] = im1
                    
    # Collapsing the double complex to a single complex
    Z = P.base_ring()
    k = Z.base_ring()
    homComplex = ProjectiveComplex(base_ring=k)
    renumberingDictionary = {}
    for (i,j) in doubleComplexObjects:
        for (a,b) in doubleComplexObjects[(i,j)]:
            Source = P.objects(i)[a]
            Target = Q.objects(j)[b]
            homs = Source.hom(Target)
            for hom_index in range(0, len(homs)):
                # Add a copy of the base field for each hom in the correct degree.
                # Store the hom as a basis, along with source and target indices.
                homComplex.addObject(j-i, GradedProjectiveModuleOverField(k,
                                                                          1,
                                                                          - Z.deg(homs[hom_index]) + Target.twist() - Source.twist(),
                                                                          name="k",
                                                                          basis={((i,a), (j,b)): homs[hom_index]}))
                # Remember where the object is stored.
                renumberingDictionary[(i,j,a,b,hom_index)] = len(homComplex.objects(j-i))-1 
    # Now add maps
    for ((i,j,a,b), (I,J,A,B)) in doubleComplexMaps.keys():
        maps = doubleComplexMaps[((i,j,a,b), (I,J,A,B))]
        for (alpha, beta) in maps.keys():
            homComplex.addMap(j-i, renumberingDictionary[(i,j,a,b,alpha)], renumberingDictionary[(I,J,A,B,beta)], maps[(alpha,beta)])
    return homComplex
            

def cone(P, Q, M):
    '''
    The cone of M: P -> Q. 
    M must define a map of chain complexes from P to Q.
    '''
    if DEBUG:
        if not checkMap(P, Q, M):
            raise TypeError("Not a chain map. Cannot make a cone.")

    D = P.directSum(Q.shift(-1))
    for place in M.keys():
        for (i,j) in M.get(place,{}):
            D.addMap(place, i, j+len(P.objects(place+1)), M[place][(i,j)])

    return D
    

def checkMap(P, Q, M):
    '''
    Check that M defines a map of chain complexes from P to Q.
    M must have the type {i: d} where i is an integer and d is a dictionary {(a,b): r} whose associated matrix defines the map from P to Q (by right multiplication).
    '''

    min_index = min(P.min_index(), Q.min_index())
    max_index = max(P.max_index(), Q.max_index())
    for i in range(min_index, max_index):
        dPi = matrix(len(P.objects(i)), len(P.objects(i+1)), P.maps(i))
        dQi = matrix(len(Q.objects(i)), len(Q.objects(i+1)), Q.maps(i))
        Mi = matrix(len(P.objects(i)), len(Q.objects(i)), M.get(i,{}))
        Mip1 = matrix(len(P.objects(i+1)), len(Q.objects(i+1)), M.get(i+1,{}))
        for k, v in (dPi*Mip1 - Mi*dQi).dict().items():
            if not P.objects(i)[k[0]].is_zero(v):
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
    def is_zero(self, r):
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
    def is_zero(self, r):
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
    def qPolynomial(self, variables = [var('q')]):
        return variables[0]^self._grade

    @cached_method    
    def invert(self, r):
        if r != 0:
            return 1/r
        else:
            raise TypeError("Not invertible.")

