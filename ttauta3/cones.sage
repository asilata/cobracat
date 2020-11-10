# Order of variables: (a,b,c,x)

a = vector((1,0,0,0))
b = vector((0,1,0,0))
c = vector((0,0,1,0))
x = vector((0,0,0,1))

ieqab = vector((0,1,-1,0,0))  # a - b > 0
ieqcb = vector((0,0,-1,1,0)) # c - b > 0
ieqcbx = vector((0,0,-1,1,-1)) # c - b - x > 0
ieqbcx = vector((0,0,1,-1,-1)) # b - c - x > 0
ieqax = vector((0,1,0,0,-1))  # a - x > 0
ieqpositive = [vector((0,1,0,0,0)), vector((0,0,1,0,0)), vector((0,0,0,1,0)), vector((0,0,0,1,0)), vector((0,0,0,0,1))]

P1 = Polyhedron(ieqs=[ieqab, ieqcb, ieqcbx, ieqax] + ieqpositive) # This is a pyramid = union of following two simplices
P11 = Polyhedron(rays = [a,c,a+b+c,2*c+x+a+b])
P12 = Polyhedron(rays = [a,c,a+c+x,2*c+x+a+b])

P2 = Polyhedron(ieqs=[ieqab, ieqcb, ieqcbx, -1*ieqax] + ieqpositive) # This is already simplicial

P3 = Polyhedron(ieqs=[ieqab, ieqcb, -1*ieqcbx, ieqax] + ieqpositive) 
# This is a union of three simplicies, probably in multiple ways, but we guess the following based on avoiding intersections.
P31 = Polyhedron(rays=[a, a+b+c, a+b+c+x, a+b+2*c+x])
P32 = Polyhedron(rays=[a, a+x, a+b+c+x, a+b+2*c+x])
P33 = Polyhedron(rays=[a, a+x, a+c+x, a+b+2*c+x])

P4 = Polyhedron(ieqs=[ieqab, -1*ieqcb, ieqbcx, ieqax] + ieqpositive) # This is already simplicial

P5 = Polyhedron(ieqs=[ieqab, -1*ieqcb, -1*ieqbcx, ieqax] + ieqpositive) # This is a pyramid = union of following two simplices
P51 = Polyhedron(rays=[a, a+b+c, a+b+x, a+b+c+x])
P52 = Polyhedron(rays=[a, a+x, a+b+x, a+b+c+x])
