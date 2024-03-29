def internalTwist(ob, deg):
        newObjects = {}
        newMaps = {}
        for i in range(ob.minIndex(), ob.maxIndex()+1):
            newObjects[i] = []
            for o in ob.objects(i):
                newObjects[i].append(o.twistBy(deg))
                newMaps[i] = ob.maps(i)
        return ProjectiveComplex(ob.basering(), newObjects, newMaps, ob.names())

# * Calculate the Harder Narasimhan Filtration
def HN(ob, stab):
    '''
    Compute the HN filtration of ob with respect to the given stability condition.
    stab is a list of object ordered by increasing phases.
    The return value is the stable pieces in the HN filtration in decrease order of phase.
    '''
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
        homDegree = sorted(set([a-c for ((a,b), (c,d)) in topHom.keys()]))
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

# Return the internal twist and the heart degree of a standard stable object, given a stable object up to twist.
# Assumption: standard stable objects end in homological degree 0.
def twistShift(stable, stab):
    it = None
    heart = None
    for i in range(stable.maxIndex(), stable.minIndex()-1, -1):
        if len(stable.objects(i)) > 0:
            it = stable.objects(i)[0].twist()
            heart = i
            break
    return (it, i)
#return internalTwist(stable,-it).shift(i)

# Return the phase of the given stable object. Returns the position if there are is no list of phases given.
def phase(stable, stab, stablePhases=None):
    (it,i) = twistShift(stable, stab)
    stdStable = internalTwist(stable,-it).shift(i)
    # Currently a hack to check up to isomorphism (checks string reps)
    index = [str(s) for s in stab].index(str(stdStable))
    if stablePhases != None:
            phase = stablePhases[index]
    else:
            phase = index
    return (it-i, phase)

# Return the mass of the given stable object, given a list of masses. Returns the q-mass if a q-parameter is specified.
def mass(ob, stab, stableMasses, q=1):
    hn = HN(ob, stab)
    m = 0
    for o in hn:
        (it,i) = twistShift(o, stab)
        level = i - it
        stdStable = internalTwist(o,-it).shift(i)
        index = [str(s) for s in stab].index(str(stdStable))
        m = m + stableMasses[index]*(q**level)
    return m

# Return the support of the given stable object
def support(ob, stab):
    hn = HN(ob, stab)
    m = 0
    support = []
    for o in hn:
        (it,i) = twistShift(o, stab)
        level = i - it
        stdStable = internalTwist(o,-it).shift(i)
        index = [str(s) for s in stab].index(str(stdStable))
        if index not in support:
                support = support + [index]
    return set(support)
