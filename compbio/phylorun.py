
import optparse, shutil, sys, os, shlex
from rasmus import util


common_ext = [
    "alignext",
    "distext",
    "treeext",
    "usertreeext",
    "outputext"]

def add_common_options(o, align=True, tree=True, dist=False):
    if align:
        o.add_option("-A", "--alignext", dest="alignext",
                     metavar="ALIGN_EXT",
                     default=".align")
    if dist:
        o.add_option("-D", "--distext", dest="distext",
                     metavar="DIST_EXT",
                     default=".dist")
    if tree:
        o.add_option("-T", "--treeext", dest="treeext",
                     metavar="TREE_EXT",
                     default=".tree")
        o.add_option("-U", "--usertreeext", dest="usertreeext",
                     metavar="USER_TREE_EXT")
        o.add_option("-b", "--boot", dest="boot",
                     metavar="BOOT_ITERS",
                     type="int",
                     default=0)
        o.add_option("--no-opttree", dest="opttree", action="store_false",
                     default=True)


    o.add_option("-e", "--extra", dest="extra", metavar="EXTRA_ARGS",
                 help="extra arguments to pass to program")
    o.add_option("-O", "--outputext", dest="outputext",
                 metavar="OUTPUT_DIR_EXT",
                 default=".output")
    o.add_option("--seqtype", dest="seqtype", metavar="SEQ_TYPE",
                 default="dna")
    o.add_option("-v", "--verbose", dest="verbose",
                 action="store_true",
                 default=False,
                 help="verbose output")
    
def parse_common_options(o):
    conf, files = o.parse_args()

    if conf.extra:
        conf.extra = shlex.split(conf.extra)
    
    return conf, files


def get_basename(filename, conf, exts=common_ext):

    for extname in exts:
        try:
            ext = getattr(conf, extname)
            if ext:
                if filename.endswith(ext):
                    return util.replace_ext(filename, ext, "")
        except:
            pass
    raise Exception("file '%s' has unknown extension" % filename)

def make_output_dir(dirname):

    if os.path.exists(dirname):
        shutil.rmtree(dirname)
    os.makedirs(dirname)
