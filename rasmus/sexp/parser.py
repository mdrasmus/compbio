
import os

try:
    from rasmus.ply import lex, yacc
except ImportError:
    from ply import lex, yacc

#=============================================================================
# classes

class Sym (str):
    """S-expression symbol"""
    def __repr__(self):
        return self

#=============================================================================
# grammar



literals = ['(',')', "[", "]", "{", "}"]
t_ignore = " \t\r\n"
t_ignore_COMMENT = r";[^\n]*\n"

tokens = (
    "SYMBOL",
    "STRING",
    "NUMBER",
    "BOOL"
)


def quote(s):
    return s.replace('\\', '\\\\').replace(r'"', r'\"')

def unquote(s):
    # fix, use proper unquote
    return s.replace(r'\"', r'"').replace('\\\\', '\\')

def t_STRING(t):
    r'"(\\\\|\\"|[^"\\]+)*"'
    t.value = unquote(t.value[1:-1])
    return t

def t_SYMBOL(t):
    r'[^()\[\]{}" \#\t\r\n0-9;][^()\[\]{}" \t\r\n]*'
    t.value = Sym(t.value)
    return t

def t_NUMBER(t):
    r"[+-]?(\d+\.?|\.\d)(\d+([eE][+-]?\d+)?)?"
    if ("." not in t.value and
        "e" not in t.value and
        "E" not in t.value):
        t.value = int(t.value)
    else:
        t.value = float(t.value)
    return t

def t_BOOL(t):
    r"\#t|\#f"
    t.value = (t.value == "#t")
    return t

def t_error(t):
    raise TypeError("Unknown text '%s'" % (t.value,))


def p_sexp(p):
    """
    sexp : "(" item_list ")"
    """
    p[0] = p[2]

def p_item_list(p):
    """
    item_list : item item_list
              |
    """
    if len(p) > 1:
        p[0] = [p[1]]
        p[0].extend(p[2])
    else:
        p[0] = []

def p_item(p):
    """
    item : sexp
         | SYMBOL
         | STRING
         | NUMBER
         | BOOL
    """
    p[0] = p[1]


def p_error(p):
    if p:
        raise Exception("Syntax error at '%s'" % p.value)
    else:
        raise Exception("Syntax error")


#lex.lex()
#yacc.yacc()

outdir = os.path.dirname(__file__)
lex.lex(debug=0, optimize=1, lextab="sexp_lex", outputdir=outdir)
yacc.yacc(debug=0, optimize=1, tabmodule="sexp_tab", outputdir=outdir)

#outdir = os.path.dirname(__file__)
#lex.lex(debug=1, optimize=0, lextab="sexp_lex", outputdir=outdir)
#yacc.yacc(debug=1, optimize=0, tabmodule="sexp_tab", outputdir=outdir)


parse = yacc.parse




if __name__ == "__main__":
    print yacc.parse('(+ 223 36.6 (aaa) bbb \"ccc\")')
    print yacc.parse(r'(if (> 2 var) (cons a b) (display "no \\ \"quoted\" \\ " 22) )')
    #print yacc.parse("(7:subject(3:ref5:alice6:mother))")
