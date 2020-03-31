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
