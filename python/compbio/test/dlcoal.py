# test dlcoal

from rasmus.common import *
from rasmus import stats
from rasmus.testing import *

from compbio import coal, dlcoal
reload(coal)
reload(dlcoal)

if 1:
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


if 1:
    draw_tree_names(coal_tree, scale=.5e-7)
    print exp(dlcoal.prob_dlcoal_recon_topology(
        coal_tree, ex["coal_recon"],
        ex["locus_tree"], ex["locus_recon"], ex["locus_events"],
        ex["daughters"], fly_stree, n, duprate, lossrate,
        pretime=None, premean=None,
        maxdoom=20, nsamples=100,
        add_spec=True))

if 1:
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
