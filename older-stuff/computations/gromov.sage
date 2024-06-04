def xy(r, theta):
    return vector((r*cos(theta), r*sin(theta)))

argP2 = pi*random()
argP1 = random()*argP2
argP3 = 0

P1 = xy(random(), argP1)
P2 = xy(random(), argP2)
P3 = xy(random(), argP3)

def threeDiagonals(P1, P2, P3):
    P12 = P1 + P2
    P123 = P1 + P2 + P3
    P23 = P2 + P3
    return (((P1).norm() + (P3).norm()).n(20), 
            ((P2).norm() + (P123).norm()).n(20), 
            ((P12).norm() + (P23).norm()).n(20))


print(threeDiagonals(P1,P2,P3))

## The following was originally written in another project but transported here.
## This is the matrix for the convex stability condition 12345
M = matrix([[1,0,0,0,1,1,0,0,1,0],
            [1,1,0,0,0,0,1,0,0,1],
            [0,1,1,0,0,1,0,1,0,0],
            [0,0,1,1,0,0,1,0,1,0],
            [0,0,0,1,1,0,0,1,0,1],
            [0,1,0,0,1,1,1,0,1,1],
            [1,0,1,0,0,1,1,1,0,1],
            [0,1,0,1,0,1,1,1,1,0],
            [0,0,1,0,1,0,1,1,1,1],
            [1,0,0,1,0,1,0,1,1,1]])

var(['p' + str(a) + str(b) for b in range(1,6) for a in range(1,b)])

objs = [p12,p23,p34,p45,p15,
        p13,p24,p35,p14,p25]

## This is the matrix for the nonconvex stability condition where you move 5 in to the triangles 123 and 134.
M2 = matrix([[1,0,0,0,1,1,0,0,1,0],
             [1,1,0,0,0,0,1,0,0,1],
             [0,1,1,0,0,1,0,1,0,0],
             [0,0,1,1,0,0,1,0,1,0],
             [0,1,0,0,1,1,1,0,1,1],
             [1,0,1,0,0,1,1,1,0,1],
             [0,1,0,1,0,1,1,1,1,0],
             [1,0,1,1,1,1,1,0,0,0],
             [1,0,0,1,0,1,0,1,1,1],
             [0,0,1,0,1,0,1,1,1,1]])
