import copy

from rasmus import util

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
def func_timer(func):
    def wrapper(*args, **kwargs):
        util.tic(func.__name__)
        result = func(* args, **kwargs)
        util.toc()
        return result
    return wrapper



def func_defaults(** defaults):

    def helper(func):
        def wrapper(*args, **kwargs):
            print "HERE", func, defaults

            for k, v in defaults.items():
                if k not in kwargs:
                    kwargs[k] = copy.copy(v)

            result = func(* args, **kwargs)
            return result
        return wrapper
    return helper



if __name__ == "__main__":
    @func_timer
    def dowork(x):
        w = [1,2,3]
        print x


    dowork(99)
    
    @func_defaults(x=[1,2,3], y=7)
    def doit(x=0, y=0):
        x.append(4)
        print x, y
    
    doit()
    doit()
    
    #help(dowork)



