###############################################################################
# Timer class for timing nested sections of code

import os, sys, traceback

from util import *


class Timer:
    def __init__(self, stream = sys.stderr, maxdepth=1e1000):
        self.reset()
        self.streams =  [(stream, maxdepth)]
        self.showErrors = True
        self.showWarnings = True
        self.quiets = 0

    def start(self, msg = ""):
        t = time.time()
        if msg != "":
            self.indent()
            self._write("BEGIN %s:\n" % msg)
        self.msg.append(msg)
        self.starts.append(t)
        self.flush()
    
    def time(self):
        return self.starts[-1] - time.clock()
    
    def stop(self):
        duration = time.time() - self.starts.pop()
        msg = self.msg.pop()
        if msg != "":
            self.indent()
            self._write("END   %s: [%.3f s]\n" % (msg, duration))
        self.flush()
        return duration
    
    def log(self, *text):
        self.indent()
        for i in text:
            self._write("%s " % str(i))
        self._write("\n")
        self.flush()        
    
    def logExact(self, text):
        self._write(text)
        self.flush()
    
    def warn(self, text, offset=0):
        filename, lineno, func, code = traceback.extract_stack()[-2-offset]
        filename = os.path.basename(filename)
        
        if self.showWarnings:
            self.indent()
            self._write("WARNING: %s, line %d: %s\n" % (filename, lineno, text))
            self.flush()
        
    def error(self, text, offset=0):
        filename, lineno, func, code = traceback.extract_stack()[-2-offset]
        filename = os.path.basename(filename)
        
        if self.showErrors:
            self.indent()
            self._write("ERROR: %s, line %d: %s\n" % (filename, lineno, text))
            self.flush()
    
    def indent(self):
        for i in range(self.depth()):
            self._write("  ")
    
    def reset(self):
        self.msg = []
        self.starts = []
    
    def depth(self):
        return len(self.msg)
    
    def _write(self, text):
        for stream, maxdepth in self.streams:
            if self.depth() < maxdepth and \
               self.quiets == 0:
                stream.write(text)
    
    def write(self, text):
        self._write(text)
        self.flush()
    
    def flush(self):
        for stream, maxdepth in self.streams:
            stream.flush()
    
    def addStream(self, stream, maxdepth=1e1000):
        self.streams.append((stream, maxdepth))
    
    def removeStream(self, stream):
        indices = findneq(stream, cget(self.streams, 0))
        self.streams = sublist(self.streams, indices)

    def suppress(self):
        self.quiets += 1
    
    def unsuppress(self):
        self.quiets = max(self.quiets - 1, 0)


def globalTimer():
    if not "timer" in GLOBALS():
        GLOBALS()["timer"] = Timer()
    return GLOBALS()["timer"]

def log(*text):
    return globalTimer().log(*text)

def logger(*text):
    return globalTimer().log(*text)
        
def logExact(text):
    return globalTimer().logExact(text)

def tic(msg = ""):
    return globalTimer().start(msg)

def toc():
    return globalTimer().stop()

def indent():
    return globalTimer().indent()

def warn(text, offset=0):
    return globalTimer().warn(text, offset+1)

def error(text, offset=0):
    return globalTimer().error(text, offset+1)



def note(*text):
    print >>notefile(), " ".join(text)

def noteflush():
    return notfile().flush()

def notefile(out = None):
    if out == None:
        out = file("/dev/null", "w")
    if "notes" not in GLOBALS():
        GLOBALS()["notes"] = out
    return GLOBALS()["notes"]



################################################################################
# debugging info
#

def current_file(offset=0, abbrv=True):
    filename, lineno, func, code = traceback.extract_stack()[-2-offset]
    if abbrv:
        filename = os.path.basename(filename)
    return filename
    
def current_line(offset=0):
    filename, lineno, func, code = traceback.extract_stack()[-2-offset]
    return lineno

def current_func(offset=0):
    filename, lineno, func, code = traceback.extract_stack()[-2-offset]
    return func

def current_code(offset=0):
    filename, lineno, func, code = traceback.extract_stack()[-2-offset]
    return code



