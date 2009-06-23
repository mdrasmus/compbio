
import sys, re
import parser
from parser import quote, Sym

    

def parse(string, env=None):
    """Parse an S-expression from a string"""

    sexp = parser.parse(string)
    if env:
        sexp = process(sexp, env)
    return sexp


def read(infile, env=None):
    """Read an S-expression from a file"""
    return parse(infile.read(), env)

def process(sexp, env):
    """Evaluate an S-expression using an environment dictionary"""

    # process sexp
    if isinstance(sexp, list):
        func = env.get(sexp[0], None)
        if func and isinstance(sexp[0], Sym):
            return func(sexp, env)
        else:
            # default recursion
            sexp = [process(elm, env) for i, elm in enumerate(sexp)]

    return sexp


def prepare(obj, env=None):
    """Prepare an object heirarchy as an S-expression"""

    if env:
        for objtype, func in env:
            if isinstance(obj, objtype):
                return func(obj, env)
        else:
            if hasattr(obj, "__iter__"):
                return [prepare(elm, env) for elm in obj]
            else:
                return obj
    else:
        return obj

def write(sexp, out=sys.stdout):
    """Write an S-expression to a file"""

    if isinstance(sexp, Sym):
        out.write(sexp)

    elif isinstance(sexp, basestring):
        out.write('"' + quote(sexp) + '"')

    elif isinstance(sexp, bool):
        if sexp:
            out.write("#t")
        else:
            out.write("#f")

    elif isinstance(sexp, (int, float)):
        out.write(repr(sexp))

    elif hasattr(sexp, "__iter__"):
        # write list

        out.write("(")
        size = len(sexp)
        it = iter(sexp)
        try:
            a = it.next()
            while True:
                write(a, out)
                a = it.next()
                out.write(" ")
        except StopIteration:
            pass
                
        out.write(")")

    else:
        raise Exception("item with unknown type in sexp '%s'" % repr(sexp))


def write_pretty(sexp, out=sys.stdout, 
                 _prefix=0, _first=True):
    """Write an S-expression to a file"""

    if _first:
        out.write(" " * _prefix)

    if isinstance(sexp, Sym):
        out.write(sexp)
        return len(sexp)

    elif isinstance(sexp, basestring):
        s = '"' + quote(sexp) + '"'
        out.write(s)
        return len(s)

    elif isinstance(sexp, bool):
        if sexp:
            out.write("#t")
        else:
            out.write("#f")
        return 2

    elif isinstance(sexp, (int, float)):
        s = repr(sexp)
        out.write(s)
        return len(s)

    elif hasattr(sexp, "__iter__"):
        # write list

        i = 0
        out.write("(")
        size = 1
        size2 = 0
        it = iter(sexp)
        _prefix += 1
        try:
            a = it.next()
            while True:
                if i < 1:
                    sexp_size = write_pretty(a, out, 
                                             _prefix=_prefix, _first=False)
                    i += 1
                    size += sexp_size
                    _prefix += sexp_size
                    a = it.next()
                    out.write(" ")
                    size += 1
                    _prefix += 1
                else:
                    size2 = write_pretty(a, out, 
                                         _prefix=_prefix, _first=(i >= 2))
                    i += 1
                    a = it.next()
                    out.write("\n")
                    
        except StopIteration:
            pass
                
        out.write(")")
        return size + size2 + 1

    else:
        raise Exception("item with unknown type in sexp '%s'" % repr(sexp))


#=============================================================================
# datatypes

def sexp2dict(sexp, env=None):
    """Parse an sexp into a dictionary"""
    
    dct = {}
    it = iter(sexp)
    it.next() # skip first word

    for item in it:
        key, value = item
        dct[process(key, env)] = process(value, env)

    return dct


def dict2sexp(dct, env=None, tag="dict"):
    
    lst = [Sym(tag)]
    for key, value in dct.iteritems():
        lst.append([prepare(key, env), prepare(value, env)])
    return lst


#=============================================================================
# testing

if __name__ == "__main__":
    from StringIO import StringIO

    print parse('(+ 223 36.6 (aaa) bbb \"ccc\")')
    s = parse(r'''(if (> 2 var) (cons a b) 
                      (display "no \"quoted\" " 22 (#f #t)) )''')
    print s

    t = StringIO()
    write(s, t)
    print t.getvalue()

    write_pretty(s)
    print

    print parse('''(hello there () (dict (a 1) ("bb" 222.0) (8 #f)
                                      (more (dict (u v))))
                    ("dict" (a b) (c d)))''')
    print parse('''(hello there (dict (a 1) ("bb" 222.0) (8 #f)
                                      (more (dict (u v))))
                    ("dict" (a b) (c d)))''',
                {"dict": sexp2dict})

    write_pretty(dict2sexp({"aaa": 111, 
                            True: (((22, "abc", "adcd"), 9999),
                                   "www",
                                   (Sym("hello"), [], ("hi", 
                                                   5555,
                                                   "---")
                                     )), 
                            78: False}))
    print

    o = parse('''

    (account (usename "raz")
             (started (date August 17 2005))
             (renewed (date October 5 2009))
             (score (+ 2 (* 3 7)))
             (score2 (quote (+ 2 (* 3 7))))
             )

    ''', {"account": sexp2dict, 
          "date": lambda x,e: tuple(x[1:]),
          "+": lambda x, e : sum(map(lambda i: process(i, e), x[1:])),
          "*": lambda x, e : reduce(lambda a,b: a*b,
                                    map(lambda i: process(i, e), x[1:])),
          "quote": lambda x, e: x[1]})

    print o

    write_pretty(prepare(o, [[dict, lambda x, e: dict2sexp(x, e, "account")],
                             [tuple, lambda x, e: [Sym("date")] + list(x)]]))
    print
    print


    #=====================================================================
    # tree
    def parse_node(sexp, env):
        name = sexp[1]
        data = {}
        children = []

        for x in sexp[2:]:
            if x[0] == "node":
                children.append(parse_node(x, env))
            else:
                data[x[0]] = x[1]
        return (name, data, children)


    o = parse('''
    ;tree comments
    (node "C" (dist .2) (boot 70) (species  "root")
      (node "A" (dist .1) (boot 100)) ; branch A
      (node "B" (dist .11) (boot 98)
         (node "Human1" (dist .01) (species "Human"))
         (node "Chimp2" (dist .03) (species "Chimp"))))
    ''')

    t = process(o, {"node": parse_node})
    from pprint import pprint
    pprint(o)
    print t
    write_pretty(o)
    print
