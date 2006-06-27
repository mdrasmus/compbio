import util


def id2genome(ensid):
    if ensid.startswith("ENSG"): return "human"
    elif ensid.startswith("ENSCAFG"): return "dog"
    elif ensid.startswith("ENSMUSG"): return "mouse"
    elif ensid.startswith("ENSRNOG"): return "rat"
    else:
        raise Exception("unknown ENSEMBL ID")


def genomeComposition(genomes, comp):
    counts = {}
    for genome in genomes:
        counts[genome] = 0
    for gene in comp:
        genome = id2genome(gene)
        if genome in genomes:
            counts[genome] += 1
    return counts

def componentCompositions(order, comps):
    compositions = util.Dict(1, 0)
    for comp in comps:
        counts = genomeComposition(order, comp)
        key = []
        for genome in order:
            key.append(counts[genome])
        compositions[tuple(key)] += 1
    return compositions.data

