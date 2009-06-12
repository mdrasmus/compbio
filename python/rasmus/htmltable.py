
from rasmus import tablelib, util

def format_link(linkfunc):
    """make a format function from a link function that takes a cell's item
       and returns the link url and text"""

    def func(item):
        url, text = linkfunc(item)
        return '<a href="%s">%s</a>' % (url, text)
    return func


class HtmlTable (object):
    """
    A class for creating an HTML table from a tablelib.Table
    """

    def __init__(self, table, title=None,
                 headers=None,
                 formats=None):

        self.table = table
        self.title = title

        # column headers
        if headers is None:
            # set headers from table
            self.headers = table.headers
            
        elif isinstance(headers, dict):
            # set some headers from dict

            self.headers = list(table.headers)
            for old_header, new_header in headers.items():
                self.headers[table.headers.index(old_header)] = new_header

        else:
            # set headers from list
            self.headers = list(headers)

        # formating
        if formats is None:
            # use no formating
            self.formats = [None] * len(self.headers)
            
        elif isinstance(formats, dict):
            # set some formats from dict
            self.formats = [None] * len(self.headers)
            for header, f in formats.items():
                self.formats[table.headers.index(header)] = f            
        else:
            # set formats from list
            self.formats = list(formats)


    def write(self, out, fullpage=False):
        """Write HTML table"""
        out = util.open_stream(out, "w")

        if fullpage:
            out.write("<html>")

        if self.title:
            out.write("<head><title>%s</title></head>\n" % self.title)

        out.write("<style>.tab { border-right: 1px solid #777; border-bottom: 1px solid #777;}</style>")

        if self.title is not None:
            out.write("<h1>%s</h1>" % self.title)

        # write headers
        out.write("<table cellspacing=0 style='border: 1px solid black;'>\n")
        out.write("<tr><td class='tab'><b>#</b></td>")
        for header in self.headers:
            out.write("<td class='tab'><b>%s</b></td>" % header)
        out.write("</tr>\n")

        # write rows
        for i, row in enumerate(self.table):
            out.write("<tr><td class='tab'>%d.</td>" % (i+1))
            for j, item in enumerate(util.mget(row, self.table.headers)):
                
                if self.formats[j] is not None:
                    # write formating
                    out.write("<td class='tab'>%s&nbsp;</td>" %
                               self.formats[j](item))
                else:
                    out.write("<td class='tab'><nobr>%s&nbsp;</nobr></td>" % str(item))
            out.write("</tr>\n")    
            
        out.write("</table>")

        if fullpage:
            out.write("</html>")

