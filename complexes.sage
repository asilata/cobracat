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
        # All actions are right actions!
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
            matrices[i] = matrix(targetDim, sourceDim, self._maps.get(i, {})).transpose() # Right actions!

        for i in range(smallest, largest-1):
            if matrices[i] * matrices[i+1] != 0:
                print "Differential squared not zero at " + str(i) + "."
                return False

        return True

    def directSum(self, Q):
        R = Q.copy()
        smallest = min(self._objects.keys())
        largest = max(self._objects.keys())

        for i in range(smallest, largest+1):
            R.objects(i).extend(self.objects(i))
        for k in range(smallest, largest):
            lk, lkplusone = len(Q.objects(k)), len(Q.objects(k+1))
            kmaps = self.maps(k)
            for (p,q) in kmaps:
                R.maps(i)[(p+lk,q+lkplusone)] = kmaps[(p,q)]
        return R.cleanup()
    
    def minimizeAt(self, place):
        # Assumption: all non-zero maps of degree 0 are isomorphisms
        k = self._basering.base_ring()

        
        # Find an object at i and an object at (i+1) with an isomorphism between them.
        def _findIso(place):
            for i in range(0, len(self.objects(place))):
                for j in range(0, len(self.objects(place+1))):
                    fij = self.maps(place).get((i,j), self._basering(0))
                    if fij != 0 and fij.degree() == 0:
                        return i,j, fij
            return None, None, None

        alreadyMinimized = False
        while not alreadyMinimized:
            source, target, alpha = _findIso(place)
            if source == None or target == None or alpha == None:
                print("Nothing left to minimize at " + str(place))
                alreadyMinimized = True
                continue

            # Change the maps from place to place+1

            def invert(alpha):
                try:
                    return 1/alpha
                except TypeError:
                    try:
                        return alpha.inverse()
                    except AttributeError(e):
                        raise e
            
            newMapsPlace = {}
            for i in range(0, len(self.objects(place))):
                for j in range(0, len(self.objects(place+1))):
                    if (i,j) == (source, target):
                        changeij = 0
                    else:
                        changeij = self.maps(place).get((i,target), 0) * invert(alpha) * self.maps(place).get((source,j), 0)
                    newMapsPlace[(i,j)] = self.maps(place).get((i,j), 0) - changeij


            # The maps from place-1 to place and place+1 to place+2 do not need to be changed substantially, apart from the indexing.
            # Now we update the maps
            for i in range(0, len(self.objects(place))):
                for j in range(0, len(self.objects(place+1))):
                    self._maps[place][(i,j)] = newMapsPlace[(i,j)]

            # At this point, our complex is a direct sum of F (source) -> F (target) and another complex
            # We simply drop the source and the target
            self._objects[place].pop(source)
            self._objects[place+1].pop(target)

            # and re-index as needed
            matrixAtPlace = matrix(self.maps(place)).transpose() # Right action!
            newMatrixAtPlace = matrixAtPlace.delete_columns([source]).delete_rows([target])
            self._maps[place] = newMatrixAtPlace.transpose().dict()

            matrixAtPlaceMinus1 = matrix(self.maps(place-1)).transpose()
            if matrixAtPlaceMinus1.nrows() > 0:
                newMatrixAtPlaceMinus1 = matrixAtPlaceMinus1.delete_rows([source])
                self._maps[place-1] = newMatrixAtPlaceMinus1.transpose().dict()

            matrixAtPlacePlus1 = matrix(self.maps(place+1)).transpose()
            if matrixAtPlacePlus1.ncols() > 0:
                newMatrixAtPlacePlus1 = matrixAtPlacePlus1.delete_columns([target])
                self._maps[place+1] = newMatrixAtPlacePlus1.transpose().dict()

        #Finally we do a cleanup
        self.cleanUp()
        return 

