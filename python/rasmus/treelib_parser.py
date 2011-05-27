
import os, re

try:
    from rasmus.ply import lex, yacc
except ImportError:
    from ply import lex, yacc


literals = ['+', ':', ';', '(',')', ","]
t_ignore = " \t\n"

tokens = (
    "NAME",
    "DATA"
)



t_NAME = r"[\w\-_\.]+([\w\-_\. ]*[\w\-_\.])?"
t_DATA = r"[^,;\(\)]+"

#t_DATA2 = r"^([^,;\(\)\[\]]+|[^,;\(\)\[\]]*\[[^\]]*\][^,;\(\)\[\]]*)$"
#t_DATA =  r"([^,;\(\)]+|[^,;\(\)]*\[[^\]]*\])"
#[^,;\(\)]*XXX\[[^\]]*\])"
#[^,;\(\)]+)"

def t_error(t):
    raise TypeError("Unknown text '%s'" % (t.value,))


#print re.match(t_DATA2, ";")
#print re.match("^"+t_DATA+"$", ":1")
#print re.match("^"+t_DATA+"$", ":1[hi]")
#sys.exit()

#=============================================================================

def p_tree(p):
    """
    tree : subtree ';'
    """
    p[0] = p[1]


def p_subtree(p):
    """subtree : "(" branch_set ")" NAME DATA
               | "(" branch_set ")" NAME
               | "(" branch_set ")" DATA
               | "(" branch_set ")"
               | NAME DATA
               | NAME
               | DATA
    """
    if len(p) == 6:
        p[0] = (p[2], "", p[4] + p[5])
    elif len(p) == 5:
        p[0] = (p[2], "", p[4])
    elif len(p) == 4:
        p[0] = (p[2], "", "")
    elif len(p) == 3:
        p[0] = ([], p[1], p[2])
    elif len(p) == 2:
        if ":" in p[1]:
            p[0] = ([], "", p[1])
        else:
            p[0] = ([], p[1], "")


def p_branch_set(p):
    """branch_set : subtree "," branch_set
                  | subtree
    """

    if len(p) == 2:
        p[0] = [p[1]]
    else:
        p[0] = [p[1]] + p[3]
    

def p_error(p):
    if p:
        raise Exception("Syntax error at '%s'" % p.value)
    else:
        raise Exception("Syntax error")


#lex.lex()
#yacc.yacc()

outdir = os.path.dirname(__file__)
lex.lex(debug=0, optimize=1, lextab="treelib_lex", outputdir=outdir)
yacc.yacc(debug=0, optimize=1, tabmodule="treelib_tab", outputdir=outdir)

if __name__ == "__main__":
    print yacc.parse("(sss:1.0,(abc:.2, hello there:.1):2.0,abcd:4.0);")
    print yacc.parse("((xax:1.0,bbx:2));")
    print yacc.parse("((aa:1.0,bb:2)x:33,(cc:4,dd:5):6);")



#=============================================================================
# OLD parsing code



"""
   Tree --> Subtree ";" | Branch ";"
   Subtree --> Leaf | Internal
   Leaf --> Name
   Internal --> "(" BranchSet ")" Name
   BranchSet --> Branch | Branch "," BranchSet
   Branch --> Subtree Length
   Name --> empty | string
   Length --> empty | ":" number
"""

'''
def t_FLOAT(t):
    r"[+-]?(\d+\.?|\.\d)(\d+([eE][+-]?\d+)?)?"
    t.value = float(t.value)
    return t
'''

'''
def _tree(p):
    """
    tree : branch ';'
    """
    p[0] = p[1]



def subtree(p):
    """subtree : "(" branch_set ")" NAME
               | "(" branch_set ")" 
               | NAME"""
    if len(p) == 5:
        p[0] = (p[2], p[4])
    elif len(p) == 4:
        p[0] = (p[2], "")
    else:
        p[0] = ([], p[1])


def _subtree(p):
    """subtree : "(" branch_set ")"
               | NAME"""
    if len(p) == 4:
        p[0] = (p[2], "")
    else:
        p[0] = ([], p[1])


def _branch_set(p):
    """branch_set : branch "," branch_set
                  | branch 
    """

    if len(p) == 2:
        p[0] = [p[1]]
    else:
        p[0] = [p[1]] + p[3]
    

def _branch(p):
    """
    branch : subtree DATA
           | subtree
    """

    if len(p) == 3:
        p[0] = p[1] + (p[2],)
    else:
        p[0] = p[1] + ("",)
'''


'''
        # literals
        lparen    = Literal("(").suppress()
        rparen    = Literal(")").suppress()
        colon     = Literal(":").suppress()
        semicolon = Literal(":").suppress()
        comma     = Literal(",").suppress()
        point     = Literal(".")
        e         = CaselessLiteral("E")


        # terminal rules
        # name = Word(alphanums + "_" + "-" + "." + "+")
        name_part  = Word(alphanums + "_" + "-" + "." + "+")
        name_part2 = Forward()        
        name_part2 << name_part + Optional(Word(" ") + name_part2)
        name = Combine(name_part + Optional(Word(" ") + name_part2))
        fnumber = Combine(Word("+-"+nums, nums) + 
                          Optional(point + Optional(Word(nums))) +
                          Optional(e + Word("+-"+nums, nums)))
        dist      = fnumber
        bootstrap = fnumber


        # recursive rules
        subtree = Forward()
        subtreelist = Forward()
        
        subtree << \
            Group(
                (
                    (lparen + subtreelist + rparen).setResultsName("subtree") |
                    name.setResultsName("name")
                ) +
                Optional(
                    CharsNotIn(",);").setResultsName("data")
                )
            )
        subtreelist << subtree + Optional(comma + subtreelist)


        # top level rule
        tree = subtree + Word(";").suppress()
'''
