vertices = [1,2,3,4]
edges = {"a":(1,2), "b":(2,3), "c":(3,4), "d":(4,1), "x": (1,3), "y": (2,4)}

supports = [tuple(s) for s in Subsets(edges.keys(),4) if not ("x" in s and "y" in s) and ("x"in s or "y" in s)]

def states(support):
    incidences = {v:[e for e in support if v in edges[e]] for v in vertices}
    unpruned =  [(support, (l1,l2,l3,l4)) for l1 in incidences[1] for l2 in incidences[2] for l3 in incidences[3] for l4 in incidences[4]]
    pruned = [(x,y) for (x,y) in unpruned if len(uniq(sorted(y))) < 4] # must have a dumb-bell 
    return pruned

def stateToString(state):
    return ''.join(state[0]) + "-" + ''.join(state[1])

[[stateToString(state) for state in states(support)] for support in supports]
