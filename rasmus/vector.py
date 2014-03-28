
from math import sqrt


def vdot(u, v):
    assert len(u) == len(v)

    tot = 0.0
    if hasattr(u, "keys"):
        for i in u:
            tot += u[i] * v[i]
    else:
        for i in xrange(len(u)):
            tot += u[i] * v[i]

    return tot


def vproj(u, v):
    return vmuls(v, vdot(u, v) / vmag(v)**2)


def vadd(u, v):
    assert len(u) == len(v)
    if hasattr(u, "keys"):
        return map(lambda i, j: u[i] + v[j], u, v)
    else:
        return map(lambda a, b: a + b, u, v)


def vsub(u, v):
    assert len(u) == len(v)
    if hasattr(u, "keys"):
        return map(lambda i, j: u[i] - v[j], u, v)
    else:
        return map(lambda a, b: a - b, u, v)


def vmul(u, v):
    assert len(v) == len(u)
    if hasattr(u, "keys"):
        return map(lambda i, j: u[i] * v[j], u, v)
    else:
        return map(lambda a, b: a * b, u, v)


def vdiv(u, v):
    assert len(v) == len(u)
    if hasattr(u, "keys"):
        return map(lambda i, j: u[i] / float(v[j]), u, v)
    else:
        return map(lambda a, b: a / float(b), u, v)


def vidiv(u, v):
    assert len(v) == len(u)
    if hasattr(u, "keys"):
        return map(lambda i, j: u[i] / v[j], u, v)
    else:
        return map(lambda a, b: a / b, u, v)


def vmag(v):
    tot = 0.0
    for i in v:
        tot += i*i
    return sqrt(tot)


def vdist(u, v):
    tot = 0.0
    for i in xrange(len(v)):
        tot += (u[i] - v[i])**2
    return sqrt(tot)


#=============================================================================
# vector and scalar

def vadds(u, s):
    if hasattr(u, "keys"):
        return [u[i] + s for i in u]
    else:
        return [a + s for a in u]


def vsubs(u, s):
    if hasattr(u, "keys"):
        return [u[i] - s for i in u]
    else:
        return [a - s for a in u]


def vmuls(u, s):
    if hasattr(u, "keys"):
        return [u[i] * s for i in u]
    else:
        return [a * s for a in u]


def vdivs(u, s):
    s = float(s)
    if hasattr(u, "keys"):
        return [u[i] / s for i in u]
    else:
        return [a / s for a in u]


def vidivs(u, s):
    if hasattr(u, "keys"):
        return [u[i] / s for i in u]
    else:
        return [a / s for a in u]


def in_left_halfspace2(a, b, p):
    """Returns True is point p is to the left of line a<->b.
       where left is defined as standing at a and facing towards b"""
    return (b[0]-a[0]) * (p[1]-a[1]) - (b[1]-a[1]) * (p[0]-a[0]) <= 0


def in_triangle2(a, b, c, pos):
    """Returns True is pos in triangle a,b,c"""

    clockwise = in_left_halfspace2(b, a, c)
    if clockwise:
        return (in_left_halfspace2(b, a, pos) and
                in_left_halfspace2(c, b, pos) and
                in_left_halfspace2(a, c, pos))
    else:
        return (in_left_halfspace2(a, b, pos) and
                in_left_halfspace2(b, c, pos) and
                in_left_halfspace2(c, a, pos))


def in_polygon2(pts, pos):
    """Returns True if point 'pos' in convex polygon with points 'pts'"""

    assert len(pts) >= 3, Exception("polygon must have at least 3 points")
    clockwise = in_left_halfspace2(pts[1], pts[0], pts[2])

    if clockwise:
        for i in xrange(1, len(pts)):
            if not in_left_halfspace2(pts[i], pts[i-1], pos):
                return False
        return in_left_halfspace2(pts[0], pts[-1], pos)
    else:
        for i in xrange(0, len(pts)-1):
            if not in_left_halfspace2(pts[i], pts[i+1], pos):
                return False
        return in_left_halfspace2(pts[-1], pts[0], pos)
