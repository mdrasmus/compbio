
import unittest

from rasmus import tablelib
from StringIO import StringIO


class Test (unittest.TestCase):

    def test_guess_types(self):

        text = """\
name	num	text	truth
john	-10	1a	true
matt	123	3b	true
alex	456	2c	true
mike	789	1d	false
"""

        tab = tablelib.read_table(StringIO(text))
        self.assertEquals(tab[0],
                          {'text': '1a',
                           'num': -10.0,
                           'name':
                           'john',
                           'truth': True})

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


'''




    #################################################
    # catch parse error
    if 0:
        text="""\
##types:str	int	int
name	num	num
matt	123	0
alex	456
mike	789	1
"""

        tab = readTable(StringIO.StringIO(text))
        tab.sort()

        print repr(tab)
        print tab.defaults
        print tab
        print tab.cget('name', 'num')


    #################################################
    # timing
    if 0:
        from rasmus import util

        text=["##types:" + "int\t" * 99 + "int",
              "\t".join(map(str, range(100))) ]

        for i in range(10000):
            text.append("1\t" * 99 + "1")
        text = "\n".join(text)

        stream = StringIO.StringIO(text)

        util.tic("read table")
        tab = readTable(stream)
        util.toc()


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
