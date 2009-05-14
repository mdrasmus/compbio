
import optparse, shutil, sys, os
from rasmus import util


common_ext = [
    "alignext",
    "treeext",
    "usertreeext",
    "outputext"]

def add_common_options(o):
    o.add_option("-A", "--alignext", dest="alignext",
                 metavar="ALIGN_EXT",
                 default=".align")
    o.add_option("-T", "--treeext", dest="treeext",
                 metavar="TREE_EXT",
                 default=".tree")
    o.add_option("-U", "--usertreeext", dest="usertreeext",
                 metavar="USER_TREE_EXT")
    o.add_option("-O", "--outputext", dest="outputext",
                 metavar="OUTPUT_DIR_EXT",
                 default=".output")
    o.add_option("-b", "--boot", dest="boot",
                 metavar="BOOT_ITERS",
                 type="int",
                 default=0)
    o.add_option("--seqtype", dest="seqtype", metavar="SEQ_TYPE",
                 default="dna")
    o.add_option("--no-opttree", dest="opttree", action="store_false",
                 default=True)
    
def parse_common_options(o):
    conf, files = o.parse_args()    
    return conf, files


def get_basename(filename, conf, exts=common_ext):

    for extname in exts:
        ext = getattr(conf, extname)
        if ext:
            if filename.endswith(ext):
                return util.replace_ext(filename, ext, "")
    raise Exception("file '%s' has unknown extension" % filename)

def make_output_dir(dirname):

    if os.path.exists(dirname):
        shutil.rmtree(dirname)
    os.makedirs(dirname)
