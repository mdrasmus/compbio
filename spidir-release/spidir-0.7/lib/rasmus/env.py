import os, copy



# paths variable
datapaths = ["."]


class PathError (Exception):
    """Error for when a file cannot be found in a path list"""
    pass


def findFile(filename, pathlist = datapaths, cwd=True):
    """Searches for a filename in a list of paths"""

    for path in pathlist:
        filename2 = os.path.join(path, filename)
        if os.path.exists(filename2):
            return filename2
    raise PathError("'%s' cannot be found in paths" % filename)


def addPaths(* paths, **keywords):
    """Add paths to the pathlist
    
       Custom pathlist can be specified by 'pathlist' keyword argument
       """
    
    if "pathlist" in keywords:
        pathlist = keywords["pathlist"]
    else:
        # use global datapaths
        pathlist = datapaths
    
    for path in paths:
        pathlist.extend(path.split(":"))


def addEnvPaths(varname, pathlist=datapaths):
    newpaths = os.getenv(varname)
    if newpaths != None:
        addPaths(newpaths, pathlist=pathlist)

