from Numeric import *
from LinearAlgebra import *

from rasmus import util, gnuplot

def hermite(x0, x1, v0, v1, t):
    a = array([1, t, t**2, t**3])
    M = array([[1, 0, 0, 0],
               [1, 1, 1, 1],
               [0, 1, 0, 0],
               [0, 1, 2, 3]])
    p = transpose([[x0, x1, v0, v1]])
    
    return dot(dot(a, inverse(M)), p)

def bezier2(x0, p1, p2, x1, t):
    a = array([1, t, t**2, t**3])
    M = array([[1, 0, 0, 0],
               [1, 1, 1, 1],
               [0, 1, 0, 0],
               [0, 1, 2, 3]])
    B = array([[1, 0, 0, 0],
               [0, 0, 0, 1],
               [-3, 3, 0, 0],
               [0, 0, -3, 3]])
    p = transpose([[x0, p1, p2, x1]])

    C = dot(inverse(M),B)
    
    return dot(dot(a, C), p)

def lerp(x0, x1, t):
    return x0 * (1-t) + x1 * t

def bezier(x0, x1, x2, x3, t):
    x00 = lerp(x0, x1, t)
    x01 = lerp(x1, x2, t)
    x02 = lerp(x2, x3, t)
    x10 = lerp(x00, x01, t)
    x11 = lerp(x01, x02, t)
    x20 = lerp(x10, x11, t)
    return x20

def catmullrom(x0, x1, x2, x3, t):
    x00 = lerp(x0, x1, t+1)
    x01 = lerp(x1, x2, t)
    x02 = lerp(x2, x3, t-1)
    x10 = lerp(x00, x01, (t+1)/2.0)
    x11 = lerp(x01, x02, t/2.0)
    x20 = lerp(x10, x11, t)
    return x20

def interpolate(points):
    pts = points[:1] + points + points[-1:]
    def func(t):
        t += 1
        it = int(t) - 1
        return catmullrom(pts[it], pts[it+1],
                          pts[it+2], pts[it+3], t - int(t))
    return func


def bbasis(i, p, t, times):
    if p == 0:
        if times[i] <= t and t < times[i+1]:
            return 1
        else:
            return 0
    else:
        return (t - times[i]) / (times[i+p] - times[i]) * \
               bbasis(i, p-1, t, times) + \
               (times[i+p+1]-t) / (times[i+p+1]-times[i+1]) * \
               bbasis(i+1, p-1, t, times)

def bspline(pts, times, t):
    p = len(times) - len(pts) - 1
    tot = 0
    for i in range(len(pts)):
        tot += pts[i] * bbasis(i, p, t, times)
    return tot

def bspline2(pts, times, t):
    p = len(times) - len(pts) - 1
    b = []
    for i in range(len(pts)):
        b.append(bbasis(i, p, t, times))
    return b


f = lambda x: hermite(0,0,1,-1,x)
g = lambda x: bezier(0,1,1,0,x)
h = lambda x: catmullrom(-3,0,0,-3,x)

if False:
    p = gnuplot.plotfunc(f,0,1,.1)
    p.plotfunc(g,0,1,.1)
    p.plotfunc(h,0,1,.1)

if False:
    pts = [1,2,2,3,3,1,-5,-4,-5,-3,-5,-2,-5,0]
    p = gnuplot.plotfunc(interpolate(pts), 0, len(pts)-1, .1)

def frange(n):
    return map(lambda x: x/float(n), range(n+1))

pts = [0, 4, 0]
pts2 = [0, 4, 6]

p = gnuplot.Gnuplot()

def make(n):
    times = frange(n)
    f = lambda p: lambda x: bspline(p, times, x)
    coords = map(f(pts), frange(10))
    coords2 = map(f(pts2), frange(10))

    p.plot(coords, coords2)

times = frange(10)
f = lambda p: lambda x: bspline2(p, times, x)
coords = map(f(pts), frange(100))

plots = util.unzip(coords)
for plot in plots:
    p.plot(plot, style="lines")

#make(4)
#make(5)
#make(6)
#make(7)
#make(8)

#p = util.plotfunc(f, 0, 1, .01)

