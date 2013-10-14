
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
    cmd = " ".join(["pep8"] + filenames)
    print cmd
    pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    lines = [line for line in pipe.stdout if key(line)]
    pipe.wait()
    return lines


def pyflakes_filter(line):
    """
    Standard filter for pyflakes.
    """

    return True


def pep8_filter(line):
    """
    Standard filter for pep8.
    """

    return True


def get_python_scripts(*paths):
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

        # Return filenames ending in *.py
        if filename.endswith(".py"):
            yield filename
            continue

        # Return filenames containing 'python' in the first line
        with open(filename) as infile:
            line = infile.readline()
            if "python" in line and "python-i" not in line:
                yield filename


_python_files = list(get_python_scripts(
    "rasmus",
    "compbio",
    "bin",
    "test",
    "test/rasmus",
    "test/compbio",
))


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
