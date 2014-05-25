
import os
import subprocess


def run(cmd, shell=True):
    """
    Run a command and check the return code.
    """
    pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, shell=shell)
    out = pipe.stdout.read()
    retcode = pipe.wait()
    if retcode != 0:
        print out
    return retcode


def run_pyflakes(filenames, key=lambda line: True):
    """
    Run pyflakes and return all errors.
    """
    cmd = " ".join(["pyflakes"] + filenames)
    print cmd
    pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    lines = [line for line in pipe.stdout if key(line)]
    pipe.wait()
    return lines


def run_pep8(filenames, key=lambda line: True):
    """
    Run pep8 and return all errors.
    """
    options = []
    ignore = []

    # E265 block comment should start with '# '
    ignore.append('E265')

    # E226 missing whitespace around arithmetic operator
    ignore.append('E226')

    options.extend(['--ignore', ','.join(ignore)])

    cmd = " ".join(["pep8"] + options + filenames)
    print cmd
    pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    lines = [line for line in pipe.stdout if key(line)]
    pipe.wait()
    return lines


def pyflakes_filter(line):
    """
    Standard filter for pyflakes.
    """
    ignore = [
        'from rasmus.timer import *',
        'from timer import *',
        'from rasmus.vector import *',
        'from vector import *',
        'from rasmus.plotting import *',
        'from plotting import *',
    ]
    for text in ignore:
        if text in line:
            return False

    return True


def pep8_filter(line):
    """
    Standard filter for pep8.
    """
    return True


def get_python_scripts(paths, exclude=[]):
    """
    Return the python scripts in a directory
    """
    filenames = []
    for path in paths:
        files = sorted(os.listdir(path))
        filenames.extend(os.path.join(path, filename) for filename in files)
    for filename in filenames:
        # Skip directories
        if not os.path.isfile(filename):
            continue

        # Skip files that start with exclude prefix
        if any(filename.startswith(prefix)
               for prefix in exclude):
            continue

        # Return filenames ending in *.py
        if filename.endswith(".py"):
            yield filename
            continue

        # Return filenames containing 'python' in the first line
        with open(filename) as infile:
            line = infile.readline()
            if "python" in line and "python-i" not in line:
                yield filename


_python_paths = [
    'rasmus',
    'compbio',
    'bin',
    'test',
    'test/rasmus',
    'test/compbio',
]
_ignore_files = [
    'bin/',
    'compbio/blast.py',
    'compbio/gff.py',
    'compbio/go.py',
    'compbio/mrbayes.py',
    'compbio/muscle.py',
    'compbio/nexus.py',
    'compbio/paml.py',
    'compbio/pfam.py',
    'compbio/phylip.py',
    'compbio/phylo.py',
    'compbio/phylogenomics.py',
    'compbio/phylorun.py',
    'compbio/phyml.py',
    'compbio/regionlib.py',
    'compbio/seqlib.py',
    'compbio/synteny/__init__.py',
    'compbio/synteny/fuzzy.py',
    'compbio/synteny/strict.py',
    'compbio/vis/argvis.py',
    'compbio/vis/transsvg.py',
    'rasmus/common.py',
    'rasmus/ply/',
    'rasmus/treelib_tab.py',
    'rasmus/treelib_lex.py',
    'rasmus/sexp/sexp_tab.py',
    'rasmus/sexp/sexp_lex.py',
    'rasmus/vis/',
]
_python_files = list(get_python_scripts(_python_paths, exclude=_ignore_files))


def test_pyflakes():
    """
    Run pyflakes on python code base.
    """
    lines = run_pyflakes(_python_files, key=pyflakes_filter)

    if len(lines) > 0:
        print "pyflakes errors:"
        print "".join(lines)
        raise Exception()


def test_pep8():
    """
    Ensure pep8 compliance on python code base.
    """
    lines = run_pep8(_python_files, key=pep8_filter)

    if len(lines) > 0:
        print "pep8 errors:"
        print "".join(lines)
        raise Exception()
