point_vars = [var(x) for x in ['s1','s2','s3','sx','t1','t2','t3','tx']]
vector_vars = [var(x) for x in ['a1','b1','a2','b2','a3','b3','ax','bx']]

s1 = -3
s2 = -2
s3 = -1
t1 = 1
t2 = 1
t3 = 2
sx = 7
tx = 10

v1 = vector([a1,b1])
v2 = vector([a2,b2])
v3 = vector([a3,b3])
vx = vector([ax,bx])
p1 = vector([s1,t1])
p2 = vector([s2,t2])
p3 = vector([s3,t3])
px = vector([sx,tx])

eqs = [
    v1.dot_product(p1),
    v2.dot_product(p2),
    v3.dot_product(p3),
    vx.dot_product(px),
    (v1+vx).dot_product(p1+px),
    (v2+vx).dot_product(p2+px),
    (v3+vx).dot_product(p3+px),
    (v1+v2+vx).dot_product(p1+p2+px),
    (v1+v3+vx).dot_product(p1+p3+px),
    (v2+v3+vx).dot_product(p2+p3+px),
    (v1+v2+v3+vx).dot_product(p1+p2+p3+px),
    (v1+v2+v3+2*vx).dot_product(p1+p2+p3+2*px)
]

def eq_to_ieq(expr):
    return [0] + [expr.coefficient(x) for x in vector_vars]

ineqs = [eq_to_ieq(e) for e in eqs]

P = Polyhedron(ieqs = ineqs, base_ring=QQ)
