
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
                 format=None):

        self.table = table
        self.title = title

        if headers is None:
            self.headers = table.headers
        else:
            self.headers = headers
        
        if format is None:
            self.format = [None] * len(self.headers)
        else:
            self.format = format            


    def write(self, out):
        """Write HTML table"""
        out = util.openStream(out, "w")

        out.write("<style>.tab { border-right: 1px solid #777; border-bottom: 1px solid #777;}</style>")

        if self.title is not None:
            out.write("<h1>%s</h1>" % self.title)

        # write headers
        out.write("<table cellspacing=0 style='border: 1px solid black;'>")    
        out.write("<tr><td class='tab'><b>#</b></td>")
        for header in self.headers:
            out.write("<td class='tab'><b>%s</b></td>" % header)
        out.write("</tr>")

        # write rows
        for i, row in enumerate(self.table):
            out.write("<tr><td class='tab'>%d.</td>" % (i+1))
            for j, item in enumerate(util.mget(row, self.table.headers)):
                
                if self.format[j] is not None:
                    # write formating
                    out.write("<td class='tab'>%s&nbsp;</td>" %
                               self.format[j](item))
                else:
                    out.write("<td class='tab'>%s&nbsp;</td>" % str(item))
            out.write("</tr>")    
            
        out.write("</table>")
