import os
import sys


import util



class Latex:
    def __init__(self, filename):
        self.filename = filename
        self.out = file(filename, "w")
        self.tempfiles = []
    
    
    def __del__(self):
        for f in self.tempfiles:
            if os.path.isfile(f):
                os.remove(f)
    
    
    def begin(self, options="", text=None):
        if text == None:
            text = r"""
\documentclass[%s]{article}
\usepackage{graphicx}  %% figures
\usepackage{multicol}  %% columns
\usepackage{latexsym}  %% gnuplot symbols

%% page formating
\pdfpagewidth 8.5in
\pdfpageheight 11in 
\textwidth=6.5in
\textheight=9in
\oddsidemargin=0in
\evensidemargin=0in
\topmargin=0in
\headheight=0in
\headsep=0in
\footskip=0.5in
\parindent=0.0cm
\parskip=0.3cm


%%\def\Diamond{\bullet}


\begin{document}
               
               
               """ % options
        
        self.out.write(text)
    
        
    def end(self):
        self.out.write("\n\n\\end{document}")
        self.out.close()
    
    
    def beginMulticols(self, cols):
        self.out.write("\n\n\\begin{multicols}{%d}\n" % cols)
    
    def endMulticols(self):
        self.out.write("\n\n\\end{multicols}\n")
    
    
    def text(self, text):
        self.out.write(text)
    
    
    def graphic(self, image, width=3):
        template = r"""
        
        %%\begin{center}
        \includegraphics[width=%fin]{%s}
        
        %%\end{center}
        
        """
        
        self.out.write(template % (width, image))
    
    
    def plot(self, plot, width=3):
        scalew = width / 5.
        scaleh = scalew
        
        #plot.gnuplot("set size %f, %f\n" % (scalew, scaleh))

        tmpfile = util.tempfile(".", "latex_", ".ps")
        self.tempfiles.append(tmpfile.replace(".ps", ".pdf"))
        plot.save(tmpfile)
        if os.path.exists(tmpfile):
            os.system("ps2pdf %s" % tmpfile)
            os.remove(tmpfile)
            print "removed", tmpfile

            self.graphic(tmpfile.replace(".ps", ".pdf"), width=width)
            self.out.write("\\vspace{-.5in}\n\n")
        else:
            util.warn("could not create '%s'" % tmpfile)
        
        
        #tmpfile = util.tempfile(".", "latex_", ".ps")
        #self.tempfile.appedn(tmpfile)
        
        
        #plot.save(tmpfile)
        #self.out.write()
        #self.out.write("\n\n")
    
    
    def makePdf(self):
        os.system("pdflatex %s" % self.filename)
    
    def view(self):
        os.system("pdflatex %s" % self.filename)
        os.system("xpdf %s" % self.filename.replace(".tex", ".pdf"))



if __name__ == "__main__":
    tex = Latex("out.tex")
    tex.begin()
    
    tex.text(r"my name is $\alpha \beta$.")
    
    tex.end()
    tex.view()
    





