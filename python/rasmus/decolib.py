import util

def decorator(deco):
    def new_decorator(f):
        g = deco(f)
        g.__name__ = f.__name__
        g.__doc__ = f.__doc__
        g.__dict__.update(f.__dict__)
        return g
    # Now a few lines needed to make simple_decorator itself
    # be a well-behaved decorator.
    new_decorator.__name__ = deco.__name__
    new_decorator.__doc__ = deco.__doc__
    new_decorator.__dict__.update(deco.__dict__)
    return new_decorator


@decorator
def timer_deco(func):
    def wrapper(*args, **kwargs):
        util.tic(func.__name__)
        result = func(* args, **kwargs)
        util.toc()
        return result
    return wrapper

if __name__ == "__main__":
    @timer_deco
    def dowork(x):
        w = [1,2,3]
        print x


    dowork(99)
    help(dowork)
