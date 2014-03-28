

import sys, os, copy
import tempfile as temporaryfile


#=============================================================================
# R plotting

# private global rplot state
_rplot_pdf = None
_rplot_temp = False
_rplot_viewer = "xpdf"





class LazyR (object):
    """Allows lazy loading of rpy"""

    def __init__(self, name):
        self.__name = name
        self.__thread = None
        self.__setup = False
        self.__rpy_version = 0
        self.__rpy = None

    def _process_events(self):
        if self.__rpy_version == 2 and self.__thread is None:
            import rpy2.rinterface as rinterface
            import time
            import threading

            def r_refresh(interval=0.05):
                while True:
                    try:
                        rinterface.process_revents()
                    except:
                        pass
                    time.sleep(interval)

            self.__thread = threading.Timer(0.1, r_refresh)
            self.__thread.start()

    def _setup(self):
        try:
            import rpy
            self.__rpy = rpy
            self.__rpy_version = 1
        except:
            import rpy2.rpy_classic as rpy
            rpy.set_default_mode(rpy.BASIC_CONVERSION)
            self.__rpy = rpy
            self.__rpy_version = 2
            self._process_events()
        globals()[self.__name] = rpy.r        

    def __getattr__(self, attr):
        if not self.__setup:
            self._setup()
        return self.__rpy.r.__getattr__(attr)

    def __call__(self, *args, **kargs):
        if not self.__setup:
            self._setup()        
        return self.__rpy.r(*args, **kargs)


rp = LazyR("rp")



def rplot_start(filename, *args, **kargs):
    """Starts a new PDF file"""
    
    global _rplot_pdf
    if "useDingbats" not in kargs:
        kargs["useDingbats"] = False
    rp.pdf(file=filename, *args, **kargs)
    _rplot_pdf = filename

def rplot_end(show=False):
    """Ends a PDF file"""
    
    global _rplot_pdf, _rplot_temp
    rp.dev_off()
    
    if show:
        if _rplot_temp:
            os.system("('%s' '%s'; rm '%s') &" % 
                      (_rplot_viewer, _rplot_pdf, _rplot_pdf))
        else:
            os.system("('%s' '%s') &" % (_rplot_viewer, _rplot_pdf))

    _rplot_pdf = None
    _rplot_temp = False


def rplot(func, *args, **kargs):
    """
    Wrapper for plotting with R.

    Makes sensible plot labels and manages pdf plotting
    """
    
    global _rplot_pdf, _rplot_temp
    
    kargs.setdefault("xlab", "")
    kargs.setdefault("ylab", "")
    kargs.setdefault("main", "")
    
    # parse my args
    if "pdf" in kargs:
        _rplot_pdf = kargs["pdf"]
        _rplot_temp = True  
        del kargs["pdf"]
        rp.pdf(file=_rplot_pdf)
        self_open = True
    else:
        self_open = False


    if "show" in kargs:
        show = kargs["show"]
        del kargs["show"]
    else:
        show = False           

    # prepare tempfile if needed
    if _rplot_pdf is None:
        f, _rplot_pdf = temporaryfile.mkstemp(".pdf", "rplot_")
        _rplot_temp = True
        os.close(f)        
        rp.pdf(file=_rplot_pdf)
        
        # force show for tempfile
        self_open = True
        show = True
            
    
    if "pdf_close" in kargs:
        close = kargs["pdf_close"]
        del kargs["pdf_close"]
    else:
        if self_open:
            close = True
        else:
            close = False
    
    
    # make R call   
    rp.__getattr__(func)(*args, **kargs)
    
    # close PDF and show
    if close:
        if _rplot_pdf is not None:
            rplot_end(show)


def rplotfunc(self, cmd, func, start, end, step, **options):
    """Plots a function using R"""
    x = []
    y = []
    
    while start < end:
        try:
            y.append(func(start))
            x.append(start)
        except ZeroDivisionError:
            pass
            start += step
    rplot(cmd, x, y, **options)
    

def rhist(*args, **kargs):
    """Plots a histogram"""
    rplot("hist", *args, **kargs)

def rplot_set_viewer(viewer):
    global _rplot_viewer
    _rplot_viewer = viewer

def rplot_get_viewer(viewer):
    return _rplot_viewer
