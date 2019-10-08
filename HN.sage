# * Calculate the Harder Narasimhan Filtration

def HN(ob, stab):
    '''
    Compute the HN filtration of ob with respect to the given stability condition.
    stab is a list of object ordered by increasing phases.
    The return value is the stable pieces in the HN filtration in decrease order of phase.
    '''
    def internalTwist(ob, deg):
        newObjects = {}
        newMaps = {}
        for i in range(ob.minIndex(), ob.maxIndex()+1):
            newObjects[i] = []
            for o in ob.objects(i):
                newObjects[i].append(o.twistBy(deg))
                newMaps[i] = ob.maps(i)
        return ProjectiveComplex(ob.basering(), newObjects, newMaps, ob.names())

    def findTop(ob):
        highestHeartFound = -Infinity
        highestObjectFound = None
        highestHomFound = None
        for stable in stab:
            H = hom(stable, ob)
            H.minimize()
            for i in range(H.minIndex(), H.maxIndex()+1):
                for j in range(0, len(H.objects(i))):
                    if H.objects(i)[j].grade()-i >= highestHeartFound: # Found a better object mapping in
                        highestHeartFound = (H.objects(i)[j].grade()-i)
                        highestObjectFound = internalTwist(stable, H.objects(i)[j].grade()).shift(-i)
                        highestHomFound = H.objects(i)[j].basis()
        return (highestObjectFound, highestHomFound)

    def isZero(ob):
        return all(len(ob.objects(i)) == 0  for i in range(ob.minIndex(), ob.maxIndex()+1))
            
    HNFiltration = []
    while not isZero(ob):
        (topOb, topHom) = findTop(ob)
        homDegree = uniq([a-c for ((a,b), (c,d)) in topHom.keys()])
        assert len(homDegree) == 1
        deg = homDegree[0]
        HNFiltration.append(topOb)
        M = {}
        for ((a,b), (c,d)) in topHom.keys():
            if c not in M.keys():
                M[c] = {}
            M[c][(b,d)] = M.get(c,{}).get((b,d), 0) + topHom[(a,b),(c,d)]
        ob = cone(topOb, ob, M).shift(1)
        ob.minimize()

    return HNFiltration
    
