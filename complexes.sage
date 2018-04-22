class ProjectiveComplex(object):
    def __init__(self, objects = {}, maps = {}):
        # Put checks to make sure the maps are well-defined
        # and they form a complex.
        self._objects = objects.copy()
        self._maps = maps.copy()
        self._names = {}

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

    def addObject(self, place, obj, name = None):
        if place not in self._objects:
            self._objects[place] = []
        self._objects[place].append(obj)

        if name != None:
            self._names[hash(obj)] = name
        
        
    def addMap(self, place, i, j, scalar):
        # Add sanity checks
        if place not in self._maps:
            self._maps[place] = {}
        self._maps[place][(i,j)] = scalar

    def checkComplexity(self):
        smallest = min(self._objects.keys())
        largest = max(self._objects.keys())

        for k in self._maps.keys():
            if k >= largest or k < smallest:
                print "Rogue map at place " + str(k) + "."
                return False

        matrices = {}
        for i in range(smallest, largest):
            sourceDim = len(self._objects.get(i, []))
            targetDim = len(self._objects.get(i+1, []))
            try:
                matrices[i] = matrix(targetDim, sourceDim, self._maps.get(i, {}))
            except IndexError(e):
                print "Index out of range in map at " + str(i) + "."
                return False

        for i in range(smallest, largest-1):
            if matrices[i+1] * matrices[i] != 0:
                print "Differential squared not zero at " + str(i) + "."
                return False

        return True
