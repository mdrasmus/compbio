
from StringIO import StringIO
from unittest import TestCase

from rasmus.sexp import Sym
from rasmus.sexp import dict2sexp
from rasmus.sexp import parse
from rasmus.sexp import prepare
from rasmus.sexp import process
from rasmus.sexp import sexp2dict
from rasmus.sexp import write
from rasmus.sexp import write_pretty


class SexpTests(TestCase):
    def test_parse(self):
        exp = parse('(+ 223 36.6 (aaa) bbb \"ccc\")')
        expected = [Sym('+'), 223, 36.6, [Sym('aaa')], Sym('bbb'), 'ccc']
        self.assertEqual(exp, expected)

        exp = parse(r'''(if (> 2 var) (cons a b)
                         (display "no \"quoted\" " 22 (#f #t)) )''')
        expected = [
            Sym('if'),
            [Sym('>'), 2, Sym('var')],
            [Sym('cons'), Sym('a'), Sym('b')],
            [Sym('display'), 'no "quoted" ', 22, [False, True]]]
        self.assertEqual(exp, expected)

        exp = parse('''(hello there () (dict (a 1) ("bb" 222.0) (8 #f)
                                          (more (dict (u v))))
                        ("dict" (a b) (c d)))''')
        expected = [Sym('hello'), Sym('there'), [], [Sym('dict'), [Sym('a'), 1], ['bb', 222.0], [8, False], [Sym('more'), [Sym('dict'), [Sym('u'), Sym('v')]]]], ['dict', [Sym('a'), Sym('b')], [Sym('c'), Sym('d')]]]  # nopep8
        self.assertEqual(exp, expected)

    def test_write(self):
        exp = parse(r'''(if (> 2 var) (cons a b)
                         (display "no \"quoted\" " 22 (#f #t)) )''')
        t = StringIO()
        write(exp, t)
        expected = '(if (> 2 var) (cons a b) (display "no \\"quoted\\" " 22 (#f #t)))'  # nopep8
        self.assertEqual(t.getvalue(), expected)

        expected = '''\
(if (> 2
       var)
    (cons a
          b)
    (display "no \\"quoted\\" "
             22
             (#f #t)))'''
        t = StringIO()
        write_pretty(exp, t)
        self.assertEqual(t.getvalue(), expected)

    def test_eval(self):
        exp = parse('''(hello there (dict (a 1) ("bb" 222.0) (8 #f)
                                          (more (dict (u v))))
                        ("dict" (a b) (c d)))''',
                    {"dict": sexp2dict})
        expected = [Sym('hello'), Sym('there'),
                    {8: False,
                     Sym('more'):
                     {Sym('u'):
                      Sym('v')},
                     'bb': 222.0,
                     Sym('a'): 1},
                    ['dict', [Sym('a'), Sym('b')], [Sym('c'), Sym('d')]]]
        self.assertEqual(exp, expected)

    def test_dict2sexp(self):
        exp = dict2sexp({"aaa": 111,
                         True: (((22, "abc", "adcd"), 9999),
                                "www",
                                (Sym("hello"), [], ("hi",
                                                    5555,
                                                    "---")
                                 )),
                         78: False})
        expected = [
            Sym('dict'),
            [True, (((22, 'abc', 'adcd'), 9999),
                    'www', (Sym('hello'), [], ('hi', 5555, '---')))],
            ['aaa', 111],
            [78, False]
        ]
        self.assertEqual(exp, expected)

    def test_eval2(self):
        o = parse('''
        (account (usename "raz")
                 (started (date August 17 2005))
                 (renewed (date October 5 2009))
                 (score (+ 2 (* 3 7)))
                 (score2 (quote (+ 2 (* 3 7))))
                 )

        ''', {"account": sexp2dict,
              "date": lambda x, e: tuple(x[1:]),
              "+": lambda x, e: sum(map(lambda i: process(i, e), x[1:])),
              "*": lambda x, e: reduce(lambda a, b: a*b,
                                       map(lambda i: process(i, e), x[1:])),
              "quote": lambda x, e: x[1]})

        expected = {Sym('score2'): [Sym('+'), 2, [Sym('*'), 3, 7]], Sym('usename'): 'raz', Sym('score'): 23, Sym('renewed'): (Sym('October'), 5, 2009), Sym('started'): (Sym('August'), 17, 2005)}  # nopep8
        self.assertEqual(o, expected)

        write_pretty(
            prepare(o, [[dict, lambda x, e: dict2sexp(x, e, "account")],
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
