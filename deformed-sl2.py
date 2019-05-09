import matplotlib.pyplot as plt
q = 1
z0 = 1
N = 100000
pointsDict = {}

def s1(z):
    if z is None:
        return 1
    elif q+z == 0:
        return None
    else:
        return z/float(q+z)

def t2(z):
    if z is None:
        return None
    else:
        return 1 + z/q

def constructPoints():
    c = 1
    while c < N:
        if c == 1:
            pointsDict[1] = z0
            c = c + 1
            continue
        z = pointsDict[c // 2]
        if c % 2 == 0:
            pointsDict[c] = s1(z)
        else:
            pointsDict[c] = t2(z)
        c = c+1

constructPoints()
vs = pointsDict.values()
ks = [0 for x in vs]
plt.plot(vs,ks, 'r.')
plt.plot([1-q],[0],'bo')
if q != 1:
    plt.plot([q/float(1-q)],[0],'bo')
plt.axis([-N.bit_length(), N.bit_length(), -0.1, 0.1])
plt.show(block=False)


