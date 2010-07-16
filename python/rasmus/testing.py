



def fequal(f1, f2, rel=.0001, eabs=1e-12):
    """assert whether two floats are approximately equal"""
    
    if f1 == f2:
        return
    
    if f2 == 0:
        err = f1
    elif f1 == 0:
        err = f2
    else:
        err = abs(f1 - f2) / abs(f2)
    x = (err < rel)

    if abs(f1 - f2) < eabs:
        return

    assert x, "%e != %e  [err=%f]" % (f1, f2, err)


