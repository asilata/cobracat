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
        return list(self._objects.get(i, []))

    def maps(self, i):
        return self._maps.get(i, {}).copy()

    def addObject(self, place, obj, name = None):
        if place not in self._objects:
            self._objects[place] = []
        self._objects[place].append(obj)

        if place not in self._maps:
            self._maps[place] = {}

        if name != None:
            self._names[hash(obj)] = name
        
    def cleanUp(self):
        for i in self._objects.keys():
            for k in self._maps[i].keys():
                if self._maps[i][k] == 0:
                    self._maps[i].pop(k)
        
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
        k = self._basering.base_ring()
        if len(self.objects(i)) == 0 or len(self.objects(i+1)) == 0:
            print("Nothing to minimize at " + str(i))
            return

        # Build a bipartite graph whose vertices are the objects at i and (i+1).
        graph = {(i,j): [(i+1,k) for k in range(0, len(self.objects(i+1)))
                         if self.maps(i).get((j,k),self._basering(0)).degree() == 0]
                 for j in range(0, len(self.objects(i))) }

        for comp in Graph(graph).connected_components():
            # Build a matrix
            print "Component: " + str(comp)
            sources = [x for x in comp if x[0] == i]
            targets = [x for x in comp if x[0] == i+1]

            if len(sources) == 0 or len(targets) == 0:
                continue

            M = matrix(k,{(b,a): self.maps(i).get((sources[a][1], targets[b][1]), 0) for a in range(0, len(sources)) for b in range(0, len(targets))})
            Minv = matrix(self._basering, M.pseudoinverse())
            print Minv
            # Change all the maps
            newMaps = {}
            for x in range(0,len(self.objects(i))):
                for y in range(0,len(self.objects(i+1))):
                    xW = matrix(self._basering, [[self.maps(i).get((x, z), 0) for (_,z) in targets]]).transpose()
                    Vy = matrix(self._basering, [[self.maps(i).get((z, y), 0) for (_,z) in sources]])
                    print xW, Vy
                    changeMap = (Vy * Minv * xW) [(0,0)]
                    oldMap = self.maps(i).get((x,y), 0)
                    newMaps[(x,y)] = oldMap - changeMap

            for x in range(0, len(self.objects(i))):
                for y in range(0, len(self.objects(i+1))):
                    self._maps[i][(x,y)] = newMaps[(x,y)]
                    
            newM = matrix(k,{(b,a): self.maps(i).get((sources[a][1], targets[b][1]), 0) for a in range(0, len(sources)) for b in range(0, len(targets))})

            # We now factor out the isomorphisms, keeping only ker(M) and ker(Minv)
            newSourceObjects = [self.objects(i)[sources[0][1]] for j in M.right_kernel().basis()]
            newTargetObjects = [self.objects(i)[targets[0][1]] for j in Minv.right_kernel().basis()]
            print newSourceObjects, newTargetObjects

