class ProjectiveObject(object):
    def __init__(self, name, twist = 0):
        self._name = name
        self._twist = twist

    def __str__(self):
        return name + "<" + str(twist) + ">"

    def __repr__(self):
        return name + "<" + str(twist) + ">"

    def twistBy(n = 1):
        self._twist = self._twist + n

class ProjectiveComplex(object):
    def __init__(self, objects = {}, maps = {}):
        # Put checks to make sure the maps are well-defined
        # and they form a complex.
        self._objects = objects
        self._maps = maps

    def addObject(self, place, obj):
        oldObjects = self._objects.get(place, [])
        oldObjects.append(obj)
        
    def addMap(self, place, i, j, scalar):
        # Add sanity checks
        self.__maps.get(place, {})[(i,j)] = scalar

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

    def minimize(self):
        #???
