
from StringIO import StringIO
import unittest

from rasmus import tablelib


class Test (unittest.TestCase):

    def test_guess_types(self):

        self.assertEquals(tablelib.guess_type('10'), int)
        self.assertEquals(tablelib.guess_type('-10'), int)
        self.assertEquals(tablelib.guess_type('10.0'), float)
        self.assertEquals(tablelib.guess_type('True'), bool)
        self.assertEquals(tablelib.guess_type('False'), bool)
        self.assertEquals(tablelib.guess_type('true'), bool)
        self.assertEquals(tablelib.guess_type('false'), bool)
        self.assertEquals(tablelib.guess_type('some text'), str)

    def test_str2bool(self):

        self.assertEquals(tablelib.str2bool('true'), True)
        self.assertEquals(tablelib.str2bool('false'), False)
        self.assertEquals(tablelib.str2bool('True'), True)
        self.assertEquals(tablelib.str2bool('False'), False)

    def test_type_names(self):

        self.assertEqual(tablelib.parse_type('int'), int)
        self.assertEqual(tablelib.parse_type('float'), float)
        self.assertEqual(tablelib.parse_type('str'), str)
        self.assertEqual(tablelib.parse_type('string'), str)
        self.assertEqual(tablelib.parse_type('bool'), bool)

        self.assertEqual(tablelib.format_type(int), 'int')
        self.assertEqual(tablelib.format_type(float), 'float')
        self.assertEqual(tablelib.format_type(str), 'string')
        self.assertEqual(tablelib.format_type(bool), 'bool')

    def test_read(self):
        text = """\
##types:str	int	float	bool
name	num	real	truth
matt	-123	10.0	true
alex	456	2.5	false
mike	789	-30.0	false
"""
        expected = [
            {'real': 10.0, 'num': -123, 'name': 'matt', 'truth': True},
            {'real': 2.5, 'num': 456, 'name': 'alex', 'truth': False},
            {'real': -30.0, 'num': 789, 'name': 'mike', 'truth': False},
        ]
        tab = tablelib.read_table(StringIO(text))
        self.assertEqual(tab, expected)

    def test_table_types(self):

        # Explict types.
        text = """\
##types:string	int	string	bool
name	num	text	truth
john	-10	1a	true
"""
        tab = tablelib.read_table(StringIO(text), guess_types=False)
        self.assertEquals(tab[0],
                          {'text': '1a',
                           'num': -10,
                           'name': 'john',
                           'truth': True})

        # Do not guess types, always use string.
        text = """\
name	num	text	truth
john	-10	1a	true
"""
        tab = tablelib.read_table(StringIO(text), guess_types=False)
        self.assertEquals(tab[0],
                          {'text': '1a',
                           'num': '-10',
                           'name': 'john',
                           'truth': 'true'})

        # Guess types from first row.
        text = """\
name	num	text	truth
john	-10	1a	true
"""
        tab = tablelib.read_table(StringIO(text))
        self.assertEquals(tab[0],
                          {'text': '1a',
                           'num': -10,
                           'name': 'john',
                           'truth': True})

        # Only specify some types with no guessing.
        tab = tablelib.read_table(StringIO(text), guess_types=False,
                                  types={'truth': bool})
        self.assertEquals(tab[0],
                          {'text': '1a',
                           'num': '-10',
                           'name': 'john',
                           'truth': True})

        # Only specify some types with guessing.
        tab = tablelib.read_table(StringIO(text),
                                  types={'truth': str})
        self.assertEquals(tab[0],
                          {'text': '1a',
                           'num': -10,
                           'name': 'john',
                           'truth': 'true'})

    def test_nheaders(self):

        text = """\
##types:str	int	int
#
# hello
#
name	0	1
matt	123	3
alex	456	2
mike	789	1
"""
        tab = tablelib.read_table(StringIO(text), nheaders=0)

        tab.add_col('extra', bool, False)
        for row in tab:
            row['extra'] = True

        self.assertEquals(set(tab[0].keys()), set([0, 1, 2, 'extra']))

    def test_sort(self):
        text = """\
##types:str	int	int
name	num	num2
matt	123	3
alex	456	0
mike	789	1
"""

        tab = tablelib.read_table(StringIO(text))
        tab.sort()
        self.assertEqual(tab.cget('name', 'num'),
                         [['alex', 'matt', 'mike'], [456, 123, 789]])

    def test_table_data(self):

        expected_tab = [
            {'a': 1, 'b': 2, 'c': 3},
            {'a': 4, 'b': 5, 'c': 6},
        ]

        expected_tab2 = [
            {0: 1, 1: 2, 2: 3},
            {0: 4, 1: 5, 2: 6},
        ]

        # parse list of dicts
        data = [
            {'a': 1, 'b': 2, 'c': 3},
            {'a': 4, 'b': 5, 'c': 6},
        ]
        tab = tablelib.Table(data)
        self.assertEqual(tab, expected_tab)

        # parse list of lists with header
        data = [
            ['a', 'b', 'c'],
            [1, 2, 3],
            [4, 5, 6],
        ]
        tab = tablelib.Table(data)
        self.assertEqual(tab, expected_tab)

        # parse list of lists with separate header
        data = [
            [1, 2, 3],
            [4, 5, 6],
        ]
        tab = tablelib.Table(data, nheaders=0, headers=['a', 'b', 'c'])
        self.assertEqual(tab, expected_tab)

        # parse list of lists with no header
        data = [
            [1, 2, 3],
            [4, 5, 6],
        ]
        tab = tablelib.Table(data, nheaders=0)
        self.assertEqual(tab, expected_tab2)

        # parse list of lists with header but no data
        data = [
            ['a', 'b', 'c'],
        ]
        tab = tablelib.Table(data)
        self.assertEqual(tab, [])

    def test_read_error(self):
        text = """\
##types:str	int	int
name	num	num
matt	123	0
alex	456	2
mike	789	1
"""
        self.assertRaises(tablelib.TableException,
                          lambda: tablelib.read_table(StringIO(text)))

        text = """\
##types:str	int	int
name	num	num
matt	123	0
alex	456	2	extra
mike	789	1
"""
        self.assertRaises(tablelib.TableException,
                          lambda: tablelib.read_table(StringIO(text)))

        text = """\
##types:str	int	int
name	num	num2
matt	123	0
alex	456	not_an_int
mike	789	1
"""
        self.assertRaises(tablelib.TableException,
                          lambda: tablelib.read_table(StringIO(text)))

        text = """\
##types:str	int	int
name	num	num2
matt	123	0
alex	456
mike	789	1
"""
        self.assertRaises(tablelib.TableException,
                          lambda: tablelib.read_table(StringIO(text)))

    def test_write(self):
        text = """\
##types:string	int	float	bool
name	num	real	truth
matt	-123	10.0	True
alex	456	2.5	False
mike	789	-30.0	False
"""

        text_no_comments = """\
name	num	real	truth
matt	-123	10.0	True
alex	456	2.5	False
mike	789	-30.0	False
"""

        # Write
        tab = tablelib.read_table(StringIO(text))
        out = StringIO()
        tab.write(out, comments=True)
        self.assertEqual(out.getvalue(), text)

        # Write no comments
        tab = tablelib.read_table(StringIO(text))
        out = StringIO()
        tab.write(out)
        self.assertEqual(out.getvalue(), text_no_comments)

        # Write string
        tab = tablelib.read_table(StringIO(text))
        out = StringIO()
        out.write(str(tab))
        self.assertEqual(out.getvalue(), text_no_comments)

        expected_repr = ("""\
name  num   real      truth  \n\
matt  -123   10.0000   True  \n\
alex   456    2.5000  False  \n\
mike   789  -30.0000  False  \n\
""")
        # Write repr
        tab = tablelib.read_table(StringIO(text))
        out = StringIO()
        out.write(repr(tab))
        self.assertEqual(out.getvalue(), expected_repr)

        # read with no header and write with header
        expected = ("""\
a	b	c
1	2	3
4	5	6
""")
        data = [
            [1, 2, 3],
            [4, 5, 6],
        ]
        tab = tablelib.Table(data, nheaders=0, headers=['a', 'b', 'c'])
        out = StringIO()
        tab.write(out, nheaders=1)
        self.assertEqual(out.getvalue(), expected)

    def test_add(self):

        # parse list of dicts
        data = [
            {'a': 1, 'b': 2, 'c': 3},
            {'a': 4, 'b': 5, 'c': 6},
            {'a': 7, 'b': 8, 'c': 9}
        ]
        tab = tablelib.Table(data)
        tab.add(a=10, b=11, c=12)
        self.assertEqual(tab[-1], {'a': 10, 'b': 11, 'c': 12})

        tab.append({'a': 13, 'b': 14, 'c': 15})
        self.assertEqual(tab[-1], {'a': 13, 'b': 14, 'c': 15})

        expected = [
            {'a': 1, 'c': 3, 'b': 2},
            {'a': 4, 'c': 6, 'b': 5},
            {'a': 7, 'c': 9, 'b': 8},
            {'a': 10, 'c': 12, 'b': 11},
            {'a': 13, 'c': 15, 'b': 14},
            {'a': 1, 'c': 3, 'b': 2},
            {'a': 4, 'c': 6, 'b': 5},
        ]
        tab.extend([
            {'a': 1, 'b': 2, 'c': 3},
            {'a': 4, 'b': 5, 'c': 6},
        ])
        self.assertEqual(tab, expected)

    def test_add_col(self):

        # parse list of dicts
        data = [
            {'a': 1, 'b': 2, 'c': 3},
            {'a': 4, 'b': 5, 'c': 6},
            {'a': 7, 'b': 8, 'c': 9}
        ]
        expected = [
            {'a': 1, 'c': 3, 'b': 2, 'd': True},
            {'a': 4, 'c': 6, 'b': 5, 'd': False},
            {'a': 7, 'c': 9, 'b': 8, 'd': False},
        ]
        tab = tablelib.Table(data)
        tab.add_col('d', bool, data=[True, False, False])
        self.assertEqual(tab, expected)

    def test_remove_col(self):

        # parse list of dicts
        data = [
            {'a': 1, 'b': 2, 'c': 3},
            {'a': 4, 'b': 5, 'c': 6},
            {'a': 7, 'b': 8, 'c': 9}
        ]
        expected = [
            {'a': 1, 'b': 2},
            {'a': 4, 'b': 5},
            {'a': 7, 'b': 8},
        ]
        tab = tablelib.Table(data)
        tab.remove_col('c')
        self.assertEqual(tab, expected)

        expected = [
            {'a': 1},
            {'a': 4},
            {'a': 7},
        ]
        tab = tablelib.Table(data)
        tab.remove_col('b', 'c')
        self.assertEqual(tab, expected)

    def test_join(self):

        tab1 = tablelib.Table([
            {'a': 1, 'b': 2, 'c': 3},
            {'a': 4, 'b': 5, 'c': 6},
            {'a': 7, 'b': 8, 'c': 9}
        ])
        tab2 = tablelib.Table([
            {'a': 1, 'd': 2, 'e': 3},
            {'a': 4, 'd': 5, 'e': 6},
            {'a': 7, 'd': 18, 'e': 19}
        ])
        tab3 = tablelib.Table([
            {'a': 1, 'b': 2, 'c': 3, 'd': 2, 'e': 3},
            {'a': 4, 'b': 5, 'c': 6, 'd': 5, 'e': 6},
            {'a': 7, 'b': 8, 'c': 9, 'd': 18, 'e': 19},
        ])
        tab4 = tablelib.Table([
            {'a': 1, 'b': 2, 'c': 3, 'd': 2, 'e': 3},
            {'a': 4, 'b': 5, 'c': 6, 'd': 5, 'e': 6},
        ])

        # join by column name
        join = tablelib.join_tables(
            (tab1, 'a', ['a', 'b', 'c']),
            (tab2, 'a', ['d', 'e']))
        self.assertEqual(join, tab3)

        # join by key function
        join = tablelib.join_tables(
            (tab1, lambda x: x['a'], ['a', 'b', 'c']),
            (tab2, lambda x: x['a'], ['d', 'e']))
        self.assertEqual(join, tab3)

        # join by key function
        join = tablelib.join_tables(
            (tab1, lambda x: (x['a'], x['b']), ['a', 'b', 'c']),
            (tab2, lambda x: (x['a'], x['d']), ['d', 'e']))
        self.assertEqual(join, tab4)

    def test_histtab(self):

        data = "aaaacbb"
        hist = tablelib.histtab(data)
        expected = tablelib.Table([
            {'item': 'a', 'count': 4, 'percent': 4 / 7.},
            {'item': 'b', 'count': 2, 'percent': 2 / 7.},
            {'item': 'c', 'count': 1, 'percent': 1 / 7.},
        ])
        self.assertEqual(hist, expected)

        data = tablelib.Table([
            ['first', 'second', 'third'],
            ['a', 'b', 1],
            ['a', 'b', 1],
            ['a', 'c', 2],
            ['c', 'a', 1],
            ['a', 'b', 2]])
        hist = tablelib.histtab(data, cols=['first', 'second'])
        expected = tablelib.Table([
            {'count': 3, 'second': 'b', 'percent': 0.6, 'first': 'a'},
            {'count': 1, 'second': 'a', 'percent': 0.2, 'first': 'c'},
            {'count': 1, 'second': 'c', 'percent': 0.2, 'first': 'a'},
        ])
        self.assertEqual(hist, expected)

'''
    #################################################
    # specialized types
    if 1:
        text="""\
##types:str	int	strand_type
name	num	strand
matt	123	+
alex	456	-
mike	789	+
john	0	+
"""
        class strand_type:
            def __init__(self, text=None):
                if text == None:
                    self.val = True
                else:
                    if text == "+":
                        self.val = True
                    elif text == "-":
                        self.val = False
                    else:
                        raise Exception("cannot parse '%s' as strand_type" %
                                        str(text))


            def __str__(self):
                if self.val:
                    return "+"
                else:
                    return "-"


        def strand_parser(text=None):
            if text == None:
                return True
            else:
                if text == "+":
                    return True
                elif text == "-":
                    return False
                else:
                    raise Exception("cannot parse '%s' as strand_type" %
                                    str(text))

        def strand_formatter(val):
            if val:
                return "+"
            else:
                return "-"

        strand_type = TableType(strand_parser, strand_formatter)


        stream = StringIO.StringIO(text)
        tab = readTable(stream, type_lookup=[["strand_type", strand_type]])
        print tab.types
        print tab

    #################################################
    # quoted strings
    if 1:
        text=\
r"""##types:str	bool	quoted_string
name	num	blah
matt	True	hello\tthere
alex	False	hello\nthere
mike	True	hello\\there
john	False	hello\n\\\nthere
"""

        stream = StringIO.StringIO(text)
        tab = readTable(stream)
        print tab.types
        print tab


    #################################################
    # python data structures/code
    if 1:
        def eval2(text=None):
            if text == None:
                return None
            else:
                return eval(text)

        python_type = TableType(eval2, str)



        tab = Table(headers=["name", "list"],
                    types={"list": python_type},
                    type_lookup=[["python", python_type]])


        tab.append({"name": "matt", "list": [1,2,3]})
        tab.append({"name": "mike", "list": [4,5,6]})
        tab.append({"name": "alex", "list": [7,8,9]})

        tab.write()

    ##################################################
    # join tables
    if 1:
        tab1 = Table([[0, 1, 2],
                      [1, 3, 4],
                      [2, 5, 6],
                      [3, 7, 8]],
                     headers=['a', 'b', 'c'])
        tab2 = Table([[0, 6, 6],
                      [1, 7, 7],
                      [3, 8, 8]],
                     headers=['a2', 'b2', 'c2'])

        tab3 = joinTables(
            (tab1, lambda x: x['a']+1, ['c', 'b']), (tab2, 'a2', ['b2']))

        print tab3
'''
