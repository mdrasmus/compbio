import scipy
from scipy.linalg import expm, logm, eig
from scipy import dot


def makeQ(bgfreq, tsvratio):
    a, c, g, t = bgfreq
    y = c + t
    r = a + g
    p = r / y
    
    b = 1. / (2*r*y * (1 + tsvratio))
    ay = (r*y*tsvratio - a*g - c*t) / \
         (2*(1.+tsvratio)*(y*a*g*p + r*c*t))
    ar = p * ay

    print "b =", b
    print "ay =", ay
    print "ar =", ar
    
    Q = scipy.array([[0, b*c, ar*g/r+b*g, b*t],
                     [b*a, 0, b*g, ay*t/y+b*t],
                     [ar*a/r+b*a, b*c, 0, b*t],
                     [b*a, ay*c/y+b*c, b*g, 0]])
    for i in xrange(4):
        tot = 0
        for j in xrange(4):
            if i != j:
                tot += Q[i][j]    
        Q[i][i] = - tot
    
    return Q


# make substitution matrix
# P = e^(Q*t)

bgfreq = [.2, .25, .3, .25]
#bgfreq = [.25, .25, .25, .25]
tsvratio = 4
time = 2
Q = makeQ(bgfreq, tsvratio)
P = expm(Q * time)

#s = dot(P, [1, 0, 0, 0])
#print s, s[1]+s[2]+s[3]

print Q
print P

x="""
0.304214 0.110134 0.280301 0.110134 
0.110134 0.304214 0.110134 0.280301 
0.420451 0.165201 0.444364 0.165201 
0.165201 0.420451 0.165201 0.444364 
"""

P2 = list2matrix(map(float, x.split()), nrows=4)

