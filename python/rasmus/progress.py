###############################################################################
# Progress classes

# TODO: add curses support?

from rasmus import util


class Progress:
    def __init__(self, *args):
       if len(args) == 1:           
           self.pos = 0
           self.end = float(args[0])
           self.step = .1
       elif len(args) == 2:
           self.pos = args[0]
           self.end = float(args[1])
           self.step = .1
       elif len(args) == 3:
           self.pos = args[0]
           self.end = float(args[1])
           self.step = args[2]
       else:
           raise Exception("wrong number of arguments")
    
       self.prog = self.pos
   
    def update(self, stream = None, msg = "progress %2.0f%%"):
        self.pos += 1
        if (self.pos > self.prog):
            self.prog += self.step * self.end
            if stream != None:
                print>>stream, msg % (100 * self.pos / self.end)
            else:
                util.log(msg % (100 * self.pos / self.end))


class ProgressBar (Progress):
    def __init__(self, *args, **dargs):
        Progress.__init__(self, *args)
        self.width = 60
        self.step = 1 / self.width
        self.bar = 0
        
        if "title" in dargs:
            title = dargs["title"]
        else:
            title = "progress"
        
        util.log("+-" + title + ("-"*(self.width-len(title)-1)) + "+")
        util.indent()
        util.logExact("|")
        self.printBar()
    
    def update(self):
        self.pos += 1        
        if (self.pos > self.prog):
            self.prog += int(self.step * self.end)
            self.printBar()
            if self.pos == self.end:
                util.logExact("|\n")
    
    def printBar(self):
        amount = int((self.pos / self.end * self.width) - self.bar)
        util.logExact("*" * amount)
        self.bar += amount


class ProgressBar (Progress):
    def __init__(self, *args, **dargs):
        Progress.__init__(self, *args)
        self.width = 60
        self.step = 1 / self.width
        self.bar = 0
        
        if "title" in dargs:
            title = dargs["title"]
        else:
            title = "progress"
        
        util.log("+-" + title + ("-"*(self.width-len(title)-1)) + "+")
        util.indent()
        util.logExact("|")
        self.printBar()
    
    def update(self):
        self.pos += 1        
        if (self.pos > self.prog):
            self.prog += int(self.step * self.end)
            self.printBar()
            if self.pos == self.end:
                util.logExact("|\n")
    
    def printBar(self):
        amount = int((self.pos / self.end * self.width) - self.bar)
        util.logExact("*" * amount)
        self.bar += amount


class FancyProgressBar (Progress):
    def __init__(self, *args, **dargs):
        Progress.__init__(self, *args)
        self.width = 60
        self.step = 1 / self.width
        self.bar = 0
        
        if "title" in dargs:
            title = dargs["title"]
        else:
            title = "progress"
        
        import curses
        self.pad = curses.initscr()
        
        self.pad.addstr(0, 0, "+-" + title + ("-"*(self.width-len(title)-1)) + "+")
        self.pad.addstr(1, 0, "|")
        self.printBar()
        self.pad.refresh()
    
    def update(self):
        self.pos += 1        
        if (self.pos > self.prog):
            self.prog += int(self.step * self.end)
            self.printBar()
            self.pad.refresh()
            
            if self.pos == self.end:
                self.pad.addstr(1, 1 + self.bar, "|\n")
                
                import curses
                curses.endwin()
    
    def printBar(self):
        amount = int((self.pos / self.end * self.width) - self.bar)
        self.pad.addstr(1, 1 + self.bar, "*" * amount)
        self.bar += amount

if __name__ == "__main__":
    import time
    
    
    util.tic("hi")
    
    prog = FancyProgressBar(100)
    for i in range(100):
        time.sleep(.01)        
        prog.update()
    util.toc()
    
    
