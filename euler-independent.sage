from itertools import combinations
n = 5
domain = range(1,n+1)

def eulerianCombinations(lst):
    paths = []
    for i in range(2, len(lst)+1):
        paths.extend([list(x) for x in combinations(lst,i)])
    return paths

def combToPath(c):
    pathDict = {}
    for i in range(c[0],c[-1]+1):
        if i == c[0] or i == c[-1]:
            pathDict[i] = 0
        elif i in c:
            pathDict[i] = 1
        else:
            pathDict[i] = -1
    return pathDict

def pathToComb(p):
    path = []
    for k in p.keys():
        if p[k] == 0 or p[k] == 1:
            path = path + [k]
    return path

eulerianPaths = [combToPath(c) for c in eulerianCombinations(domain)]

def isCrossing(p1,p2):
    commonKeys = [k for k in p1.keys() if k in p2.keys()]
    return not(all([p1[k] >= p2[k] for k in commonKeys]) or all([p1[k] <= p2[k] for k in commonKeys]))

def isNonPointed(p1,p2):
    k1,k2 = p1.keys(), p2.keys()
    if k1[0] == k2[-1] or k1[-1] == k2[0]:
        return True
    return False

def crossingOrNPPaths(p):
    for q in eulerianPaths:
        if isCrossing(p,q):
            print pathToComb(q), "crossing"
        if isNonPointed(p,q):
            print pathToComb(q), "non-pointed"

def isGood(p,qs):
    return all([not isNonPointed(p,q) and not isCrossing(p,q) for q in qs])

def maximalStates(ps):
    states = [set([])]
    maxStates = []
    while len(states) > 0:
        s = states.pop()
        spath = [combToPath(p) for p in s]
        goodPaths = [q for q in ps if q not in spath and isGood(q, spath)]
        if len(goodPaths) == 0:
            if s not in maxStates:
                maxStates.append(s)
        for q in goodPaths:
            states.append(s | set([tuple(pathToComb(q))]))
        states.sort(key=lambda x: len(x), reverse=True)
    return maxStates
    
def maxStates(ps):
    maxState = []
    maxState = maxState + pathToComb(eulerianPaths[0])
    for q in eulerianPaths:
        if isGood(q,maxState):
            maxState = maxState + [q]
    maxStates = [maxState]
    openStates = [maxState.remove(x) for x in maxState]
    closedStates = []
