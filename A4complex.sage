variables = list(var('m1','m2','m3','m4','m12','m123','m1234','m23','m234','m34'))
positivity = [v >= 0 for v in variables]

def triangle(x,y,z):
    return [x+y >= z, x+z >= y, y+z >=x]

# Assumption inequality is of the form linear-expression >= 0
def ieq_to_vector(ieq, variables):
    coefficient_list = [QQ((ieq.lhs()-ieq.rhs()).coefficient(v)) for v in variables]
    return [0]+coefficient_list

def eq_to_vector(eq, variables):
    coefficient_list = [QQ((eq.lhs()-eq.rhs()).coefficient(v)) for v in variables]
    return [0]+coefficient_list

triangles = triangle(m1,m2,m12) + triangle(m2,m3,m23) + triangle(m3,m4,m34) + triangle(m4,m1234,m123) + triangle(m1234,m1,m234) + triangle(m123,m34,m12+m4) + triangle(m12,m23,m2+m123) + triangle(m23,m34,m3+m234) + triangle(m123,m234,m1234+m23) + triangle(m12,m234,m1+m34)

ieqs = [ieq_to_vector(e, variables) for e in positivity+triangles]
P = Polyhedron(ieqs=ieqs)

def CM(m1,m2,m3,mz,mx,my):
    CMMatrix = [[0,1,1,1,1],
                [1,0,m1^2,mx^2,mz^2],
                [1,m1^2,0,m2^2,my^2],
                [1,mx^2,m2^2,0,m3^2],
                [1,mz^2,my^2,m3^2,0]]
    return matrix(CMMatrix).det()

CMConditions = [CM(m1,m2,m3,m123,m12,m23),
                CM(m2,m3,m4,m234,m23,m34),
                CM(m3,m4,m1234,m12,m34,m123),
                CM(m4,m1234,m1,m23,m123,m234),
                CM(m1234,m1,m2,m34,m234,m12)]

def random_point(s):
    return sum([randint(1,10)*r.vector() for r in s.rays()], zero_vector(s.ambient_dim()))

skel3 = []
for face3 in P.faces(4):
    center = random_point(face3.as_polyhedron())
    if all([condition.substitute({variables[i]:center[i] for i in range(0,10)}).expand() == 0 for condition in CMConditions]):
        skel3.append(face3)

def parse(eqn):
    return sum([variables[i]*eqn.vector()[i+1] for i in range(0,6)])

dualgraph = [(i,j) for i in range(0,12) for j in range(0,12) if skel2[i].as_polyhedron().intersection(skel2[j].as_polyhedron()).dim() == 2]
interior = [0,1,4,5,6,8,9,11]
skel2int = [skel2[i] for i in interior]
dualgraphint = [(i,j) for i in interior for j in interior if skel2[i].as_polyhedron().intersection(skel2[j].as_polyhedron()).dim() == 2]

# (0) P1, X, P3 
# (1) P1, P3, Z
# (2) P1 = 0, P2, P3  --- 0, 1, 4, 5, 8
# (3) P3 = 0, P1, P2
# (4) P1, P2, P3 --- 0, 2, 3, 5, 10
# (5) P1-Y-P3
# (6) P2, X, Z
# (7) Z = 0, P2, P3
# (8) P2, P1, Z
# (9) P2, P3, Z
# (10) P2 = 0, P1, P3
# (11) P2, Y, Z

# Non-convex polyhedron

trianglesNC = triangle(m1,m2,mx) + triangle(m1+m2,m3,mz) + triangle(m1,my,mz) + triangle(m2,my,m3) + triangle(m2+my,mx,mz) + triangle(m1+my,mx,m3)
ieqsNC = [ieq_to_vector(e, variables) for e in positivity+trianglesNC]
Q = Polyhedron(ieqs=ieqsNC)

skel2NC = []
for face2 in Q.faces(3):
    center = random_point(face2.as_polyhedron())
    if CMPoly.substitute({variables[i]:center[i] for i in range(0,6)}).expand() == 0:
        skel2NC.append(face2)

dualgraphNC = [(i,j) for i in range(0,9) for j in range(0,9) if skel2NC[i].as_polyhedron().intersection(skel2NC[j].as_polyhedron()).dim() == 2]        

# (0)=P2, X', Y 
# (1)=(4) P1, P2, P3 --- 0, 2, 3, 5, 10
# (2)=(10) P2 = 0, P1, P3
# (3)=Y, P1, X'
# (4)=(2) P1 = 0, P2, P3  --- 0, 1, 4, 5, 8
# (5) Y = 0, P1, P2
# (6)=(5) P1, Y, P3
# (7)=(8) P2, P1, Z
# (8)=(11) P2, Y, Z
