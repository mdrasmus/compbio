###############################################################################
# Argument parsing

from util import *   

import sys
import getopt


class OptionError (Exception):
    def __init__(self, msg):
        Exception.__init__(self)
        self.msg = msg
    def __str__(self):
        return str(self.msg)

class Option:
    def __init__(self, option_list):
        if isinstance(option_list, str):
            self.comment = option_list
            self.name = None
        elif isinstance(option_list, list):
            self.short = option_list[0]
            self.long = option_list[1]
            self.name = option_list[2]
            self.arg = option_list[3]
            self.comment = None

            # defaults
            self.flag = ":" not in self.short
            self.req = False
            self.single = False
            self.help = ""
            self.defaultGiven = False
            self.parser = lambda x: x

            # default option specifications
            if len(option_list) >= 5:
                for key, value in option_list[4].items():
                    setattr(self, key, value)

            # set default
            if "default" in dir(self):
                self.defaultGiven = True

            # set flag defaults
            if self.flag and self.single:
                self.defaultGiven = True
                self.default = False

            if self.arg.startswith("AUTO"):
                self.arg = self.arg[4:]


class Configuration (dict):
    def __init__(self, options):
        dict.__init__(self)
        self.options = options
    
    def __getitem__(self, key):
        try:
            return dict.__getitem__(self, key)
        except KeyError:
            for option in self.options:
                if option.name == key:
                    raise OptionError("Required option '%s' not given" % option.long)
            raise OptionError("Required option '%s' not given" % key)



def parseOptions(argv, options, quit=False, resthelp = ""):
    try:
        return parseArgs(argv, options, quit=quit, resthelp=resthelp, 
                         returnRest=False)
    except OptionError, e:
        print >>sys.stdout, "%s: %s" % (os.path.basename(argv[0]), e)
        sys.exit(1)



def parseArgs(argv, options, quit=False, resthelp = "", returnRest=True,
              helpOption=True):
    
    
    # add help options
    if helpOption:
        options.append(["h", "help", "help", "",
                        {"single": True,
                         "help": "display program usage"}])
    
    # setup options
    options = map(lambda x: Option(x), options)
    options2 = filter(lambda x: x.comment == None, options)
    
    # error, if no options are given
    if len(argv) < 2:
        usage(os.path.basename(argv[0]), options, resthelp)
        raise OptionError("no options given")
    
    
    # parse options
    try:
        short_opts = "".join(map(lambda x: x.short, options2))
        long_opts  = map(lambda x: x.long, options2)
        args, rest = getopt.getopt(argv[1:], short_opts, long_opts)
    except getopt.GetoptError, msg:
        usage(os.path.basename(argv[0]), options, resthelp)
        raise OptionError(msg)
    
    
    # organize options    
    lookup = {}
    for option in options2:
        lookup["-" + option.short.replace(":", "")] = option
        lookup["--" + option.long.replace("=", "")] = option
    
    
    # parse options
    conf = Configuration(options)
    given = {}
    for name, value in args:
        # figure out which option is given
        option = lookup[name]
        given[name] = 1
        
        # parse option value
        value = option.parser(value)
        
        # save option to list
        conf.setdefault(option.name, []).append(value)
    
    
    # parse options that are single or flags
    for name in given:
        option = lookup[name]
        if option.single:
            if option.flag:
                conf[option.name] = True
            else:
                conf[option.name] = conf[option.name][-1]
    
    
    # check for help option
    if helpOption and "help" in conf:
        usage(os.path.basename(argv[0]), options, resthelp)
        sys.exit(1)
        
    
    # check options
    for option in options2:
        # set default arguments
        if option.defaultGiven and option.name not in conf:
            conf[option.name] = option.default
        
        # check for required arguments
        if option.req and option.name not in conf:
            raise OptionError("required argument -%s, --%s not given" % 
                (option.short.replace(":",""), option.long.replace("=", "")))
    
    conf[""] = rest
    
    
    if returnRest:
        return conf, rest
    else:
        return conf


def usage(progname, options, resthelp = ""):
    """ 
    options = [("a:", "a_long=", "a_name", "[-a <arg>] [--a_long=<arg>]"), ...]
    """
    
    print >>sys.stderr, "Usage: %s [OPTION] %s\n" % (progname, resthelp)
    for option in options:
        if option.comment != None:
            print >>sys.stderr, option.comment
        else:
            print >>sys.stderr, "  -%s,--%s %s" % (option.short.replace(":", ""), 
                                                  option.long.replace("=", ""),
                                                  option.arg)
            if option.help != "":
                if option.help.startswith(" ") or \
                   option.help.startswith("\t"):
                    print >>sys.stderr, "%s\n" % option.help
                else:
                    print >>sys.stderr, "    %s\n" % option.help
            else:
                print >>sys.stderr
    print >>sys.stderr
    
    
"""
def getopt(* lst):
    import getopt
    param = Dict(1, [])
    
    options, rest = getopt.getopt(* lst)
    
    for option in options:
        param[option[0]].append(option[1])
    return (param.data, rest)
"""
