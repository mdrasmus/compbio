"""
 Environment variables, and paths are automatically setup with this module.
 Paths and variables differ depending on the host.

"""

import os, copy
from socket import getfqdn


hostname = getfqdn()


# paths
datapath = "/afs/csail.mit.edu/u/r/rasmus/data/"
data = datapath
genomepath = "/afs/csail.mit.edu/u/r/rasmus/data/genomes/"



# paths variable
datapaths = ["."]

class PathError (Exception):
    """Error for when a file cannot be found in a path list"""
    pass


def findFile(filename, paths = datapaths, cwd=True):
    """Searches for a filename in a list of paths"""

    for path in paths:
        filename2 = os.path.join(path, filename)
        if os.path.exists(filename2):
            return filename2
    raise PathError("'%s' cannot be found in paths" % filename)


def addPaths(path, paths=datapaths):
    """Add paths to the path variable"""
    
    datapaths[:] = path.split(":") + datapaths


def addEnvPaths(varname, paths=datapaths):
    newpaths = os.getenv(varname)
    if newpaths != None:
        addPaths(newpaths, paths)


