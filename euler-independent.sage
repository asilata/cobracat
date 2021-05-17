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
            path.append(k)
    return path

def pathToTuple(p):
    path = []
    for k in p.keys():
        if p[k] == 0 or p[k] == 1:
            path.append(k)
    return tuple(path)

trivialPaths = [domain, [domain[0],domain[-1]]]
eulerianPaths = [combToPath(c) for c in eulerianCombinations(domain) if c not in trivialPaths]
eulerianTuples = [tuple(c) for c in eulerianCombinations(domain)]

def isCrossing(p1,p2):
    commonKeys = [k for k in p1.keys() if k in p2.keys()]
    return not(all([p1[k] >= p2[k] for k in commonKeys]) or all([p1[k] <= p2[k] for k in commonKeys]))

def isNonPointed(p1,p2):
    k1,k2 = list(p1.keys()), list(p2.keys())
    if k1[0] == k2[-1] or k1[-1] == k2[0]:
        return True
    return False

compatibilityMatrix = [[0 for i in range(len(eulerianPaths))] for j in range(len(eulerianPaths))]
for i in range(len(eulerianPaths)):
    for j in range(len(eulerianPaths)):
        p,q = eulerianPaths[i], eulerianPaths[j]
        if p == q:
            compatibilityMatrix[i][j] = 0
        elif isCrossing(p,q) or isNonPointed(p,q):
            compatibilityMatrix[i][j] = 0
        else:
            compatibilityMatrix[i][j] = 1
            
def maximalStates():
    maxState = []
    maxState.append(0)
    for i in range(0,len(eulerianPaths)):
        if all(compatibilityMatrix[i][j] == 1 for j in maxState):
            maxState.append(i)

    maxStates = [sorted(maxState)]
    openEdges = [sorted(list(x)) for x in combinations(maxState, len(maxState) -1)]
    closedEdges = []

    while len(openEdges) > 0:
        e = openEdges.pop()
        closedEdges.append(e)
        print(len(openEdges), len(closedEdges))

        for i in range(0,len(eulerianPaths)):
            if all(compatibilityMatrix[i][j] == 1 for j in e):
                newMaxState = sorted(e + [i])
                if newMaxState not in maxStates:
                    maxStates.append(newMaxState)
                    newEdges = [sorted(list(x) + [i]) for x in combinations(e, len(e) -1)]
                    for n in newEdges:
                        if n not in closedEdges:
                            openEdges.append(n)
    return maxStates


