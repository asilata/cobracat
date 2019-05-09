import matplotlib.pyplot as plt
q = 1 # Value of the deformation parameter
z0 = 1 # Compute orbit of this point
N = 100000 # Apply up to N combinations of s1 and t2
pointsDict = {} # Initially empty dictionary that will hold the values of the computations.

def s1(z):
    """Computes the action of the q-matrix $\sigma_1$. Returns None if the output is the point at infinity."""
    if z is None:
        return 1
    elif q+z == 0:
        return None
    else:
        return z/float(q+z)

def t2(z):
    """Computes the action of the q-matrix $\sigma_2^{-1}. Returns None if the outpute is the point at infinity."""
    if z is None:
        return None
    else:
        return 1 + z/q

def constructPoints():
    """Populate pointsDict with (a,v) where v is the a'th computed point in the orbit. The structure is a binary tree where the left child of a node is the value of s1 and the right child of a node is the value of t2."""
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

# Run everything: populate the points, and then plot the values.        
constructPoints()
vs = pointsDict.values()
ks = [0 for x in vs]
plt.plot(vs,ks, 'r.')
plt.plot([1-q],[0],'bo')
if q != 1:
    plt.plot([q/float(1-q)],[0],'bo')
plt.axis([-N.bit_length(), N.bit_length(), -0.1, 0.1])
plt.show(block=False)


