



# d/dx f(g(x)) = f'(g(x)) g'(x)


scalar0 = ('scalar', 0.0)
scalar1 = ('scalar', 1.0)


def derivate(expr, dvar, nth=1):
    """
    Take the n'th derivative of an expression 'expr' with respect to 'dvar'
    """

    if nth == 0:
        return expr
    elif nth > 1:
        return derivate(derivate(expr, dvar, nth-1), dvar)
    elif nth < 0:
        raise Exception("nth must be >= 0")
    
    
    if expr[0] == 'scalar':
        # d/dx c = 0
        return scalar0

    elif expr[0] == 'add':
        # d/dx [f(x) + g(x)] = d/dx f(x) + d/dx g(x)
        return ('add', derivate(expr[1], dvar),
                       derivate(expr[2], dvar))

    elif expr[0] == 'mult':
        # d/dx [f(x) * g(x)] = f'(x) g(x) + f(x) g'(x)
        return ('add', ('mult', derivate(expr[1], dvar), expr[2]),
                       ('mult', expr[1], derivate(expr[2], dvar)))

    elif expr[0] == 'power' and expr[2][0] == 'scalar':
        # power rule
        # d/dx f(x)^a = a * f(x)^(a-1) * f'(x)
        a = expr[2][1]
        return ('mult', ('mult', ('scalar', a),
                                 ('power', expr[1], ('scalar', a - 1))),
                        derivate(expr[1], dvar))

    elif expr[0] == 'var':
        if expr[1] == dvar:
            # d/dx x = 1
            return scalar1
        else:
            # d/dx y = 0
            return scalar0

    # unknown: expression
    raise Exception("cannot differentiate")


def simplify(expr):

    if expr[0] == 'add':
        expr1 = simplify(expr[1])
        expr2 = simplify(expr[2])

        if expr1 == scalar0:
            return expr2
        elif expr2 == scalar0:
            return expr1
        elif expr1[0] == 'scalar' and expr2[0] == 'scalar':
            return ('scalar', expr1[1] + expr2[1])
        else:
            return ('add', expr1, expr2)

    elif expr[0] == 'mult':
        expr1 = simplify(expr[1])
        expr2 = simplify(expr[2])
        
        if expr1 == scalar1:
            return expr2
        elif expr2 == scalar1:
            return expr1
        elif expr1 == scalar0 or expr2 == scalar0:
            return scalar0
        elif expr1[0] == 'scalar' and expr2[0] == 'scalar':
            return ('scalar', expr1[1] * expr2[1])
        else:
            return ('mult', expr1, expr2)

    elif expr[0] == 'power':
        expr1 = simplify(expr[1])
        expr2 = simplify(expr[2])

        if expr2 == scalar0:
            return scalar1
        elif expr2 == scalar1:
            return expr1
        elif expr1[0] == 'scalar' and expr2[0] == 'scalar':
            return ('scalar', expr1[1] ** expr2[1])
        else:
            return ('power', expr1, expr2)

    return expr
        

def assign_vars(expr, vars):
    """Assign values to the given variables"""

    if expr[0] in ("add", "mult", "power"):
        return (expr[0], assign_vars(expr[1], vars),
                         assign_vars(expr[2], vars))
    
    elif expr[0] == "scalar":
        return expr
    
    elif expr[0] == "var":
        if expr[1] in vars:
            return ('scalar', vars[expr[1]])
        else:
            expr

    raise Expression("unkown expression")



if __name__ == "__main__":
    
    print 0, derivate(('scalar', 4), 'x')

    print "4", derivate(('mult', ('scalar', 4),
                                 ('var', 'x')), 'x')
    
    print "4", simplify(derivate(('mult', ('scalar', 4),
                                          ('var', 'x')), 'x'))

    print "4 + 2x", derivate(
        ('add',
         ('mult', ('scalar', 4), ('var', 'x')),
         ('mult', ('var', 'x'), ('var', 'x'))),
        'x')

    print "4 + 2x", simplify(derivate(
        ('add',
         ('mult', ('scalar', 4), ('var', 'x')),
         ('mult', ('var', 'x'), ('var', 'x'))),
        'x'))


    print "6x^5", derivate(
        ('power', ('var', 'x'), ('scalar', 6)),
        'x')

    print "6x^5", simplify(derivate(
        ('power', ('var', 'x'), ('scalar', 6)),
        'x'))

    
    print "4 + 6x^5", simplify(derivate(
        ('add',
         ('mult', ('scalar', 4), ('var', 'x')),
         ('power', ('var', 'x'), ('scalar', 6))),
        'x'))


    print 4 + 6 * (10**5), simplify(assign_vars(derivate(
        ('add',
         ('mult', ('scalar', 4), ('var', 'x')),
         ('power', ('var', 'x'), ('scalar', 6))),
        'x'), {'x': 10}))
