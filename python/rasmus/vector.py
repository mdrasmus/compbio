
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

def vadd(u, v):
    assert len(u) == len(v)
    if hasattr(u, "keys"):
        return map(lambda i,j: u[i] + v[j], u, v)
    else:
        return map(lambda a,b: a + b, u, v)

def vsub(u, v):
    assert len(u) == len(v)
    if hasattr(u, "keys"):
        return map(lambda i,j: u[i] - v[j], u, v)    
    else:
        return map(lambda a,b: a - b, u, v)
    
def vmul(u, v):
    assert len(v) == len(u)
    if hasattr(u, "keys"):
        return map(lambda i,j: u[i] * v[j], u, v)
    else:
        return map(lambda a,b: a * b, u, v)
    
def vdiv(u, v):
    assert len(v) == len(u)
    if hasattr(u, "keys"):
        return map(lambda i,j: u[i] / float(v[j]), u, v)    
    else:
        return map(lambda a,b: a / float(b), u, v)

def vidiv(u, v):
    assert len(v) == len(u)
    if hasattr(u, "keys"):    
        return map(lambda i,j: u[i] / v[j], u, v)
    else:
        return map(lambda a,b: a / b, u, v)


#
# vector and scalar
#

# TODO: finish using hasattr

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



