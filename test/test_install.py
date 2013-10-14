"""

Tests for installing compbio.

"""

import os

from rasmus.testing import make_clean_dir


def run_cmd(cmd):
    assert os.system(cmd) == 0


def test_install():
    """
    Test installing compbio.
    """

    make_clean_dir("test/tmp/install")
    run_cmd("python setup.py clean > /dev/null")
    run_cmd("python setup.py install --prefix=test/tmp/install > /dev/null")
    assert os.path.exists("test/tmp/install/bin/viewtree")
