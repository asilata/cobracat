# The central charges of the roots are (si,ti) for i = 1, 2, 3, x.
point_vars = [var(x) for x in ['s1','s2','s3','sx','t1','t2','t3','tx']]
# The velocity vectors of the roots are (ai,bi) for i = 1, 2, 3, x.
vector_vars = [var(x) for x in ['a1','b1','a2','b2','a3','b3','ax','bx']]
# Dummy root variables
root_vars = [var(x) for x in ['r1', 'r2', 'r3', 'rx']]

# Initialisation of the central charge.
s1,t1 = -3,1
s2,t2 = -2,1
s3,t3 = -1,1
sx,tx = 7,10

p1 = vector([s1,t1])
p2 = vector([s2,t2])
p3 = vector([s3,t3])
px = vector([sx,tx])
v1 = vector([a1,b1])
v2 = vector([a2,b2])
v3 = vector([a3,b3])
vx = vector([ax,bx])

points = [p1,p2,p3,px]
vectors = [v1,v2,v3,vx]

# Abstract list of positive roots
roots = [
    r1, r2, r3, rx,
    r1+rx, r2+rx, r3 + rx,
    r1+r2+rx, r1+r3+rx, r2+r3+rx,
    r1+r2+r3+rx, r1+r2+r3+2*rx
]

# Auxiliary list of root coefficients of positive roots in terms of the simple roots
root_coeffs = [[r.coefficient(v) for v in root_vars] for r in roots]

def coeffs_to_eq(c):
    """Convert a root coefficient to the dot product of the corresponding root with the corresponding vector."""
    v = 0
    p = 0
    for i in range(0,len(root_vars)):
        v = v + c[i]*vectors[i]
        p = p + c[i]*points[i]
    return v.dot_product(p)

# Each equation is of the form <a,v(a)> \geq 0, where a is a positive root and v(a) the velocity vector of that root.
eqs = [coeffs_to_eq(c) for c in root_coeffs]

def eq_to_ieq(expr):
    return [0] + [expr.coefficient(x) for x in vector_vars]

ieqs = [eq_to_ieq(e) for e in eqs]

P = Polyhedron(ieqs = ieqs, base_ring=QQ)

def get_state(ray):
    state = []
    for i in range(0,len(eqs)):
        eq = vector([eqs[i].coefficient(x) for x in vector_vars])
        if eq.dot_product(vector(ray)) == 0:
            state = state + [roots[i]]
    return state

def plot_central_charges():
    pts = []
    for c in root_coeffs:
        p = 0
        for i in range(0, len(root_vars)):
            p = p + c[i]*points[i]
        pts = pts + [p]

    arrows = [arrow((0,0),p, arrowsize=3, color="black") for p in pts]
    labels = [text(str(roots[i]),pts[i] + vector([0,1]), rotation="vertical", color="blue") for i in range(0,len(pts))]
    return show(sum(arrows)+sum(labels),axes_labels=['$x$','$y$'])

realStates = [[r1,  r2 + rx,  r3 + rx,  r1 + r2 + rx,  r1 + r3 + rx,  r2 + r3 + rx],
              [r2 + rx,  r3 + rx,  r1 + r2 + rx,  r1 + r3 + rx,  r2 + r3 + rx,  r1 + r2 + r3 + 2*rx],
              [r2,  r1 + rx,  r2 + rx,  r3 + rx,  r1 + r2 + rx,  r2 + r3 + rx],
              [r1 + rx,  r2 + rx,  r3 + rx,  r1 + r2 + rx,  r2 + r3 + rx,  r1 + r2 + r3 + 2*rx],
              [r3,  r1 + rx,  r2 + rx,  r3 + rx,  r1 + r3 + rx,  r2 + r3 + rx],
              [r1 + rx,  r2 + rx,  r3 + rx,  r1 + r3 + rx,  r2 + r3 + rx,  r1 + r2 + r3 + 2*rx],
              [r2,  r1 + rx,  r3 + rx,  r1 + r2 + rx,  r1 + r3 + rx,  r2 + r3 + rx],
              [r1 + rx,  r3 + rx,  r1 + r2 + rx,  r1 + r3 + rx,  r2 + r3 + rx,  r1 + r2 + r3 + 2*rx],
              [r1,  r1 + rx,  r2 + rx,  r3 + rx,  r1 + r2 + rx,  r1 + r3 + rx],
              [r1 + rx,  r2 + rx,  r3 + rx,  r1 + r2 + rx,  r1 + r3 + rx,  r1 + r2 + r3 + 2*rx],
              [r3,  r1 + rx,  r2 + rx,  r1 + r2 + rx,  r1 + r3 + rx,  r2 + r3 + rx],
              [r1 + rx,  r2 + rx,  r1 + r2 + rx,  r1 + r3 + rx,  r2 + r3 + rx,  r1 + r2 + r3 + 2*rx],
              [r1,  r3,  r2 + rx,  r1 + r2 + rx,  r1 + r3 + rx,  r2 + r3 + rx],
              [r1,  r3,  r1 + r2 + rx,  r1 + r3 + rx,  r2 + r3 + rx,  r1 + r2 + r3 + rx],
              [r1, r2, r3, r2 + rx, r1 + r2 + rx, r2 + r3 + rx],
              [r1, r2, r3, r1 + r2 + rx, r2 + r3 + rx, r1 + r2 + r3 + rx],
              [r1, r2, r3, r2 + rx, r3 + rx, r2 + r3 + rx],
              [r1, r2, r3, rx, r2 + rx, r3 + rx],
              [r1, r3, r1 + rx, r2 + rx, r3 + rx, r1 + r3 + rx],
              [r1, r3, rx, r1 + rx, r2 + rx, r3 + rx],
              [r1,  r2,  r3 + rx,  r1 + r2 + rx,  r1 + r3 + rx,  r2 + r3 + rx],
              [r1,  r2, r1 + r2 + rx,  r1 + r3 + rx,  r2 + r3 + rx,  r1 + r2 + r3 + rx],
              [r1, r3, r2 + rx, r3 + rx, r1 + r3 + rx, r2 + r3 + rx],
              [r1, r2, r2 + rx, r3 + rx, r1 + r2 + rx, r2 + r3 + rx],
              [r1, r2, r3, r3 + rx, r1 + r3 + rx, r2 + r3 + rx],
              [r1, r2, r3, r1 + r3 + rx, r2 + r3 + rx, r1 + r2 + r3 + rx],
              [r2, r3, r1 + rx, r3 + rx, r1 + r3 + rx, r2 + r3 + rx],
              [r2,  r3,  r1 + rx,  r1 + r2 + rx,  r1 + r3 + rx,  r2 + r3 + rx],
              [r2,  r3,  r1 + r2 + rx,  r1 + r3 + rx,  r2 + r3 + rx,  r1 + r2 + r3 + rx],
              [r1, r2, r3, r1 + rx, r1 + r2 + rx, r1 + r3 + rx],
              [r1, r2, r3, r1 + r2 + rx, r1 + r3 + rx, r1 + r2 + r3 + rx],
              [r1, r2, r1 + rx, r3 + rx, r1 + r2 + rx, r1 + r3 + rx],
              [r1, r2, r3, rx, r1 + rx, r3 + rx],
              [r1, r2, r3, r1 + rx, r3 + rx, r1 + r3 + rx],
              [r1, r2, r3, rx, r1 + rx, r2 + rx],
              [r1, r2, r3, r1 + rx, r2 + rx, r1 + r2 + rx],
              [r1, r3, r1 + rx, r2 + rx, r1 + r2 + rx, r1 + r3 + rx],
              [r2, r3, rx, r1 + rx, r2 + rx, r3 + rx],
              [r2, r3, r1 + rx, r2 + rx, r3 + rx, r2 + r3 + rx],
              [r2, r3, r1 + rx, r2 + rx, r1 + r2 + rx, r2 + r3 + rx],
              [r1, r2, rx, r1 + rx, r2 + rx, r3 + rx],
              [r1, r2, r1 + rx, r2 + rx, r3 + rx, r1 + r2 + rx]]

stateSets = [Set(x) for x in realStates]

def simplicesOfDim(i, states):
    answer = []
    for state in states:
        answer = answer + list(state.subsets(i))
    return list(Set(answer))

allStates = []
for i in range(1,7):
    allStates = allStates + simplicesOfDim(i, stateSets)

def link(state, states):
    return ([x for x in states if Set(state).issubset(Set(x))])

def boundary(state,states):
    return ([x for x in states if Set(x).issubset(Set(state))])

#for i in range(1,7):
#    print("Dimension " + str(i) + ": " + str(len(simplicesOfDim(i,stateSets))))

    
          
