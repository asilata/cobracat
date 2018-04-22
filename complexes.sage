class ProjectiveComplex(object):
    def __init__(self, basering, objects = {}, maps = {}):
        # Put checks to make sure the maps are well-defined
        # and they form a complex.
        self._objects = objects.copy()
        self._maps = maps.copy()
        self._names = {}
        self._basering = basering

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
                s = s + "+".join([self._names[hash(x)] if hash(x) in self._names else str(x) for x in objects])
            if i < largest:
                s = s + " → "
        return s

    def __repr__(self):
        return str(self)

    def objects(self, i):
        return self._objects.get(i, [])

    def maps(self, i):
        return self._maps.get(i, {})

    def addObject(self, place, obj, name = None):
        if place not in self._objects:
            self._objects[place] = []
        self._objects[place].append(obj)

        if name != None:
            self._names[hash(obj)] = name
        
        
    def addMap(self, place, i, j, scalar):
        if i < 0 or i >= len(self.objects(place)):
            raise IndexError("Index out of bounds")
        if j < 0 or j >= len(self.objects(place+1)):
            raise IndexError("Index out of bounds")
        
        if place not in self._maps:
            self._maps[place] = {}
        self._maps[place][(i,j)] = self._basering(scalar)

    def checkComplexity(self):
        smallest = min(self._objects.keys())
        largest = max(self._objects.keys())

        matrices = {}
        for i in range(smallest, largest):
            sourceDim = len(self._objects.get(i, []))
            targetDim = len(self._objects.get(i+1, []))
            matrices[i] = matrix(targetDim, sourceDim, self._maps.get(i, {}))

        for i in range(smallest, largest-1):
            if matrices[i+1] * matrices[i] != 0:
                print "Differential squared not zero at " + str(i) + "."
                return False

        return True
    
    def minimizeAt(self, i):
        # Assumption: all non-zero maps of degree 0 are isomorphisms
        
        if len(self.objects(i)) == 0 or len(self.objects(i+1)) == 0:
            print("Nothing to minimize at " + str(i))
            return

        sources = self.objects(i)
        targets = self.objects(i+1)

        # Build a bipartite graph whose vertices are the objects at i and (i+1).
        graph = {(i,j): [(i+1,k) for j in range(0, len(self.objects(i))) for k in range(0, len(self.objects(i+1)))
                         if self.maps(i).get((j,k),0).degree() == 0]}
        
        for comp in Graph(graph).connected_components():
            # Build a matrix
            
