# test dlcoal

import unittest

from rasmus.common import *
from rasmus import stats
from rasmus.testing import *

from compbio import coal, dlcoal
reload(coal)
reload(dlcoal)


#=============================================================================
# test prob topology from locus process


def test_dlcoal_tree(locus_tree, n, daughters, nsamples):
    """test multicoal_tree"""
    tops = {}
    
    for i in xrange(nsamples):

        # use rejection sampling
        while True:
            coal_tree, coal_recon = coal.sample_multicoal_tree(
                locus_tree, n, namefunc=lambda x: x)
            lineages = coal.count_lineages_per_branch(
                coal_tree, coal_recon, locus_tree)
            for daughter in daughters:
                if lineages[daughter][1] != 1:
                    break
            else:
                break

        
        top = phylo.hash_tree(coal_tree)
        tops.setdefault(top, [0, coal_tree, coal_recon])[0] += 1
    
    tab = Table(headers=["top", "simple_top", "percent", "prob"])
    for top, (num, tree, recon) in tops.items():
        tree2 = tree.copy()
        treelib.remove_single_children(tree2)
        tab.add(top=top,
                simple_top=phylo.hash_tree(tree2),
                percent=num/float(nsamples),
                prob=exp(dlcoal.prob_locus_coal_recon_topology(
            tree, recon, locus_tree, n, daughters)))
    tab.sort(col="prob", reverse=True)

    return tab, tops



class DLCoal (unittest.TestCase):

    def test(self):

        # params
        stree = treelib.parse_newick("((A:1000, B:1000):500, (C:700, D:700):800);")
        n = 500
        duprate = .000012
        lossrate = .000011

        # sample a locus tree with duplications
        while True:
            coal_tree, ex = dlcoal.sample_dlcoal(
                stree, n, duprate, lossrate)

            stop = False
            for d in ex["daughters"]:
                if len(d.leaves()) > 1:
                    stop = True
            if stop:
                break

        locus_tree = ex["locus_tree"]
        daughters = ex["daughters"]
        nsamples = 100
        print
        draw_tree_names(locus_tree, maxlen=8)
        
        tab, tops = test_dlcoal_tree(locus_tree, n, daughters, nsamples)
        print repr(tab[:20].get(cols=["simple_top", "percent", "prob"]))

        



#=============================================================================

if 0:
    fly_stree = parse_newick("""
(
  (
    (
      (
        (
          (
            dmel:5.32,
            (
              dsec:1.89,
              dsim:1.89
            ):3.43
          ):5.91,
          (
            dere:8.57,
            dyak:8.57
          ):2.66
        ):42.17,
        dana:53.40
      ):2.40,
      (
        dpse:1.37,
        dper:1.37
      ):54.43
    ):6.69,
    dwil:62.49
  ):1.02,
  (
    (
      dmoj:32.74,
      dvir:32.74
    ):4.37,
    dgri:37.11
  ):26.40
);""")


    def gene2species(gene):
        return "_".join(gene.split("_")[:-1])
    gen_per_myr = 1e6 / .1
    for node in fly_stree:
        node.dist *= gen_per_myr # (convert to generations)
    n = int(100e6) * 2
    duprate = .0012 / gen_per_myr
    lossrate = .0011 / gen_per_myr

    coal_tree, ex = dlcoal.sample_dlcoal(fly_stree, n, duprate, lossrate)
    #draw_tree_names(coal_tree, scale=1e-7)


if 0:
    draw_tree_names(coal_tree, scale=.5e-7)
    print exp(dlcoal.prob_dlcoal_recon_topology(
        coal_tree, ex["coal_recon"],
        ex["locus_tree"], ex["locus_recon"], ex["locus_events"],
        ex["daughters"], fly_stree, n, duprate, lossrate,
        pretime=None, premean=None,
        maxdoom=20, nsamples=100,
        add_spec=True))

if 0:
    draw_tree_names(coal_tree, scale=.5e-7, minlen=8)
    draw_tree_names(ex["locus_tree"], scale=.5e-7, minlen=8)
    print exp(dlcoal.prob_multicoal_recon_topology(
        coal_tree, ex["coal_recon"], ex["locus_tree"], n, ex["daughters"]))
    print exp(dlcoal.prob_multicoal_recon_topology2(
        coal_tree, ex["coal_recon"], ex["locus_tree"], n, ex["daughters"]))


if 0:
    draw_tree_names(coal_tree, scale=1e-7)
    print exp(coal.prob_multicoal_recon_topology(
        coal_tree, ex["coal_recon"], ex["locus_tree"], n))
    print exp(coal.prob_multicoal_recon_topology_old(
        coal_tree, ex["coal_recon"], ex["locus_tree"], n))

if 0:
    pd(coal.count_lineages_per_branch(
        coal_tree, ex["coal_recon"], ex["locus_tree"]))




#=============================================================================

show_plots = False
def show_plot():
    if show_plots:
        raw_input()


if __name__ == "__main__":

    if "--" in sys.argv:
        args = sys.argv[sys.argv.index("--")+1:]
        if "plot" in args:
            show_plots = True
        sys.argv = sys.argv[:sys.argv.index("--")]
    
    unittest.main()





