

# how long are split tracks?
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 100000     # length of locus
    k = 100        # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen

    track = []
    for j in xrange(1):
        arg = arglib.sample_arg(k, n, rho, 0, l)
        print "ARG sim"

        blocks = arglib.iter_recomb_blocks(arg)
        split_regions = defaultdict(lambda: [])

        # find regions for splits
        for block, tree in izip(blocks, arglib.iter_marginal_trees(arg)):
            for node in tree:
                if len(node.children) != 2 or node == tree.root: continue
                split = tuple(sorted(tree.leaf_names(node)))
                if len(split) == 1: continue
                split_regions[split].append(block)

        # union regions
        for split, regions in split_regions.iteritems():
            split_regions[split] =  [
                (min(r[0] for r in g), max(r[1] for r in g))
                for g in arglib.groupby_overlaps(regions)]

        for split, regions in split_regions.iteritems():
            for region in regions:
                track.append(region[1] - region[0])

    plothist(track,50)
    print "mean track", mean(track)


# visualize split tracks
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 10000      # length of locus
    k = 100        # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen
    mu = 2.5e-8    # mut/site/gen
    u = mu * l     # mut/locus/gen

    arg = arglib.sample_arg_smc(k, n, rho, 0, l)
    arg.set_ancestral()
    mut = arglib.sample_mutations(arg, u)
    mut.sort(key=lambda x: x[2])
    print "sim", len(mut)

    blocks = arglib.iter_recomb_blocks(arg)
    split_regions = defaultdict(lambda: [])
    split_nodes = defaultdict(lambda: [])

    # find regions for splits
    for block, tree in izip(blocks, arglib.iter_marginal_trees(arg)):
        for node in tree:
            if len(node.children) != 2 or node.children[0] == node.children[1]:
                continue
            split = tuple(sorted(tree.leaf_names(node)))
            if len(split) == 1: continue
            split_regions[split].append(block)
            split_nodes[split].append((block[0], node))

    # sort splits
    tracks = split_regions.items()
    tracks.sort(key=lambda x: (x[1][0][0], x[1][-1][1])) # sort by first appearance
    split_lookup = list2lookup(cget(tracks, 0))

    split_mut = set()
    mut_splits = []
    for pos, split in arglib.iter_mutation_splits(arg, mut):
        split_mut.add(split)
        mut_splits.append((pos, split))


    def click(i):
        def func():
            split = tracks[i][0]
            print i, split, split_nodes[tracks[i][0]][0]
            print [
                (min(r[0] for r in g), max(r[1] for r in g))
                for g in arglib.groupby_overlaps(tracks[i][1])]
        return func

    import summon
    #from summon.core import *
    win = summon.Window()

    # draw splits
    for i, (split, regions) in enumerate(tracks):
        col = (0, 1, 0) if split in split_mut else (1, 1, 1, .5)
        win.add_group(group(color(*col),
                group(* (lines(r[0], i, r[1], i) for r in regions)),
                hotspot("click", regions[0][0], i-.5, regions[-1][1], i+.5,
                        click(i))))

    # draw mutations
    for node, parent, pos, t in mut:
        split = tuple(sorted(x.name for x in
                             arglib.get_marginal_leaves(arg, node, pos)))
        if len(split) == 1: continue
        win.add_group(arglib.draw_mark(pos, split_lookup[split], col=(0,0,1)))

    print len(split_mut), len(split_lookup)


    # infer splits from mut_splits
    tracks2 = []
    cur = {}
    for pos, split in mut_splits:
        for split2 in cur.keys():
            if not arglib.is_split_compatible(split, split2):
                tracks2.append(((cur[split2], pos), split2))
                del cur[split2]
        if split not in cur:
            cur[split] = pos
    for split2 in cur:
        tracks2.append(((cur[split2], l), split2))

    # draw infered splits
    #for r, split in tracks2:
    #    col = (1, 0, 0)
    #    y = split_lookup[split]
    #    win.add_group(lines(color(*col), r[0], y+.2, r[1], y+.2))


# how many mutations per split track?
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 100000      # length of locus
    k = 100       # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen
    mu = 2.5e-8    # mut/site/gen
    u = mu * l     # mut/locus/gen

    nmuts = []
    for j in xrange(1):
        arg = arglib.sample_arg(k, n, rho, 0, l)
        muts = arglib.sample_mutations(arg, u)
        print "ARG sim"

        split_muts = defaultdict(lambda: 0)

        for node, parent, pos, t in muts:
            split = tuple(sorted(x.name for x in
                                 arglib.get_marginal_leaves(arg, node, pos)))
            if len(split) == 1: continue
            split_muts[split] += 1

    #rp.hist(split_muts.values(), main="", xlab="# mutations")
    plotdistrib(split_muts.values(), low=0, width=1)
    print mean(split_muts.values())


# visualize conflict graph (remove redundant conflicts)
if 0:
    rho = 1.5e-8   # recomb/site/gen
    mu = 2.5e-8    # mut/site/gen
    l = 100000     # length of locus
    k = 100        # number of lineages
    n = 2*10000    # effective popsize

    nmuts = []
    for j in xrange(1):
        arg = arglib.sample_arg(k, n, rho, 0, l)
        muts = arglib.sample_arg_mutations(arg, mu)
        print "ARG sim"

        leaves = tuple(sorted(arg.leaf_names()))
        splits = []
        for pos, split in arglib.iter_mutation_splits(arg, muts):
            splits.append(split)

        conflicts = {}

        for split in splits:
            conflicts[split] = sets.UnionFind([split])
        for split1 in splits:
            for split2 in splits:
                if not arglib.is_split_compatible_unpolar(split1, split2, leaves):
                    conflicts[split1].union(conflicts[split2])
        components = set(s.root() for s in conflicts.itervalues())
        print [len(x) for x in components]

        cmat = []
        for split1 in splits:
            cmat.append([])
            for split2 in splits:
                cmat[-1].append(int(not arglib.is_split_compatible(split1,
                                                                   split2)))

        #from rasmus import cluto
        #partids = cluto.cluster(cmat, 10, "vcluster")
        #tree = cluto.clusterTree(cmat, 10, prog="scluster")
        #perm = cluto.reorderTree(tree, cmat, prog="scluster")
        #perm = sortindex(map(sum,cmat))

        #from summon import matrix
        #v = matrix.show_dmat(cmat, rperm=perm, cperm=perm)
        v = matrix.show_dmat(cmat)


def get_unique_conflicts(splits):
    """Ignore redundant conflicts"""

    n = len(splits)
    conflicts = []
    right = [n] * n # nearest conflict to right
    left = [-1] * n  # nearest conflict to left

    leaves = tuple(sorted(arg.leaf_names()))

    # visit conflict edges in increase itervals
    for k in range(1, n):
        for i in range(n - k):
            j = i + k
            if right[i] < left[j]:
                # redundant, skip
                continue
            if not arglib.is_split_compatible_unpolar(splits[i], splits[j], leaves):
                # conflict
                conflicts.append((i, j))
                right[i] = min(right[i], j)
                left[j] = max(left[j], i)

    return conflicts, left, right


# visualize conflict graph
if 0:
    rho = 1.5e-8   # recomb/site/gen
    mu = 2.5e-8    # mut/site/gen
    l = 100000     # length of locus
    k = 100        # number of lineages
    n = 2*10000    # effective popsize


    nmuts = []
    for j in xrange(1):
        arg = arglib.sample_arg(k, n, rho, 0, l)
        muts = arglib.sample_arg_mutations(arg, mu)
        print "ARG sim"

        splits = []
        for pos, split in arglib.iter_mutation_splits(arg, muts):
            splits.append(split)

        conflicts, left, right = get_unique_conflicts(splits)
        imat = []
        for i, j in conflicts:
            imat.append((i, j, 1))
            imat.append((j, i, 1))

        from summon import matrix
        v = matrix.show_imat(imat)


# how far does one need to go along an alignment to find an incompatible mut?
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 50000      # length of locus
    k = 1000         # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen
    mu = 2.5e-8    # mut/site/gen
    u = mu * l     # mut/locus/gen

    nmuts = []
    for j in xrange(1):
        arg = arglib.sample_arg(k, n, rho, 0, l)
        muts = arglib.sample_mutations(arg, u)
        muts.sort(key=lambda x: x[2])

        blocks = arglib.iter_recomb_blocks(arg)
        splits = set()
        incompats = []

        for node, parent, pos, t in muts:
            if node.is_leaf(): continue
            split = tuple(sorted(x.name for x in
                                 arglib.get_marginal_leaves(arg, node, pos)))

            for split2 in splits:
                if not arglib.is_split_compatible(split, split2):
                    incompats.append(pos)
                    splits.clear()
                    break
            splits.add(split)

    print incompats
    lens = []
    for i in xrange(len(incompats) - 1):
        lens.append(incompats[i+1] - incompats[i])
    plothist(lens)
    print mean(lens)


#=============================================================================
# visualization

# visualize arg and mutations
if 0:
    l = 10000      # length of locus
    k = 100         # number of lineages
    n = 2*10000    # effective popsize
    rho = 1.5e-8   # recomb/site/gen
    r = rho * l    # recomb/locus/gen
    mu = 2.5e-8    # mut/site/gen
    u = mu * l     # mut/locus/gen

    arg = arglib.sample_arg(k, n, rho, 0, l)
    mut = arglib.sample_mutations(arg, u)

    def minlog(x, low=10):
        return log(max(x, low))
    #ymap = lambda x: x
    ymap = minlog

    import summon
    from summon.core import *
    win = summon.Window()
    x = 0
    step = 2
    treewidth = ilen(arg.leaves()) + step

    def trans_camera(win, x, y):
        v = win.get_visible()
        win.set_visible(v[0]+x, v[1]+y, v[2]+x, v[3]+y, "exact")

    win.set_binding(input_key("]"), lambda : trans_camera(win, treewidth, 0))
    win.set_binding(input_key("["), lambda : trans_camera(win, -treewidth, 0))

    blocks = list(arglib.iter_recomb_blocks(arg))

    for tree, block in izip(arglib.iter_marginal_trees(arg), blocks):
        pos = block[0]
        print pos

        layout = arglib.layout_arg(tree, yfunc=ymap)
        win.add_group(translate(x, 0, color(1,1,1),
                                arglib.draw_tree(tree, layout)))

        # mark responsible recomb node
        for node in tree:
            if pos != 0.0 and node.pos == pos:
                nx, ny = layout[node]
                win.add_group(arglib.draw_mark(x + nx, ny))

        # draw mut
        for node, parent, mpos, t in mut:
            if (node.name in tree and node.name != tree.root.name and
                block[0] < mpos < block[1]):
                nx, ny = layout[tree[node.name]]
                win.add_group(arglib.draw_mark(x + nx, ymap(t), col=(0,0,1)))

                for leaf in arg.leaves(node):
                    nx, ny = layout[tree[leaf.name]]
                    win.add_group(arglib.draw_mark(x + nx, ymap(0),col=(0,1,0)))

        x += treewidth

    # home
    win.home("exact")


# visualize arg and mutations
if 0:
    l = 10000      # length of locus
    k = 100         # number of lineages
    n = 2*10000    # effective popsize
    rho = 1.5e-8   # recomb/site/gen
    r = rho * l    # recomb/locus/gen
    mu = 2.5e-8    # mut/site/gen
    u = mu * l     # mut/locus/gen

    #times, events = arglib.sample_coal_recomb_times(k, n, r)
    #arg = arglib.make_arg_from_times(k, times, events)
    #arg.set_recomb_pos(0, l)
    #arg.set_ancestral()
    #arg.prune()
    arg = arglib.sample_arg(k, n, rho, 0, l)
    mut = arglib.sample_mutations(arg, u)

    #reload(arglib)
    #layout = arglib.layout_arg(arg)
    #win = arglib.show_arg(arg, mut=mut)

    cmap = rainbowColorMap(low=0, high=l)

    def minlog(x, low=10):
        return log(max(x, low))
    #ymap = lambda x: x
    ymap = minlog

    import summon
    from summon.core import *
    win = summon.Window()
    x = 0
    step = 2
    treewidth = ilen(arg.leaves()) + step

    def trans_camera(win, x, y):
        v = win.get_visible()
        win.set_visible(v[0]+x, v[1]+y, v[2]+x, v[3]+y, "exact")

    win.set_binding(input_key("]"), lambda : trans_camera(win, treewidth, 0))
    win.set_binding(input_key("["), lambda : trans_camera(win, -treewidth, 0))

    blocks = arglib.iter_recomb_blocks(arg)

    for tree, block in izip(arglib.iter_marginal_trees(arg), blocks):
        pos = block[0]
        print pos

        leaves = sorted((x for x in tree.leaves()), key=lambda x: x.name)
        layout = arglib.layout_arg(tree, leaves, yfunc=ymap)
        win.add_group(translate(x, 0, color(1,1,1,.2),
                                arglib.draw_tree(tree, layout)))

        # mark responsible recomb node
        for node in tree:
            if pos != 0.0 and node.pos == pos:
                nx, ny = layout[node]
                win.add_group(arglib.draw_mark(x + nx, ny))

        # draw mut
        for node, parent, mpos, t in mut:
            if (node.name in tree and node.name != tree.root.name and
                block[0] < mpos < block[1]):
                nx, ny = layout[tree[node.name]]
                win.add_group(arglib.draw_mark(x + nx, ymap(t), col=(0,0,1)))
            if node.name in tree and tree[node.name].parents:
                nx, ny = layout[tree[node.name]]
                py = layout[tree[node.name].parents[0]][1]
                start = arg[node.name].data["ancestral"][0][0]
                win.add_group(lines(color(0,1,0), #*cmap.get(start)),
                                    x+nx, ny, x+nx, py,
                                    color(1,1,1)))


        x += treewidth

    win.set_visible(* win.get_root().get_bounding() + ("exact",))


#=============================================================================
# understand recombinations and compatiability


# visualize compatiability
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 100000      # length of locus
    k = 60         # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen
    mu = 2.5e-8    # mut/site/gen
    u = mu * l     # mut/locus/gen

    arg = arglib.sample_arg(k, n, rho, 0, l)
    mut = arglib.sample_mutations(arg, u)
    mut.sort(key=lambda x: x[2])
    print "sim", len(mut)

    # remove uninformative muts (singletons)
    splits, mut = zip(* (a for a in ((arglib.get_mutation_split(arg, m), m)
                                  for m in mut)
                         if len(a[0]) > 1))

    # get recomb points
    rpos = arglib.get_recomb_pos(arg)
    part = []
    i = 0
    for node, parent, pos, t in mut:
        while i < len(rpos) and pos > rpos[i]:
            i += 1
        part.append(i)

    splits = [arglib.get_mutation_split(arg, m) for m in mut]

    compat = []
    for i in xrange(len(splits)):
        print "mut", i
        for j in xrange(i+1, len(mut)):
            rel = arglib.split_relation(splits[i], splits[j])
            if rel == "disjoint":
                compat.append((i, j, 1))
            elif rel in ("parent", "child", "equal"):
                compat.append((i, j, 2))
            elif rel == "conflict":
                pass
            else:
                raise Exception("unknown rel" + rel)

    import summon.matrix
    v = summon.matrix.show_imat(compat, rpart=part, cpart=part, style="quads")


# visualize redundant conflicts
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 100000      # length of locus
    k = 60         # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen
    mu = 2.5e-8    # mut/site/gen
    u = mu * l     # mut/locus/gen

    arg = arglib.sample_arg(k, n, rho, 0, l)
    mut = arglib.sample_mutations(arg, u)
    mut.sort(key=lambda x: x[2])
    print "sim", len(mut)

    # remove uninformative muts (singletons)
    splits, mut = zip(* (a for a in ((arglib.get_mutation_split(arg, m), m)
                                  for m in mut)
                         if len(a[0]) > 1))

    # get recomb points
    rpos = arglib.get_recomb_pos(arg)
    part = []
    i = 0
    for node, parent, pos, t in mut:
        while i < len(rpos) and pos > rpos[i]:
            i += 1
        part.append(i)

    splits = [arglib.get_mutation_split(arg, m) for m in mut]

    compat = []
    for i in xrange(len(splits)):
        print "mut", i
        compat.append([0] * (i+1))
        for j in xrange(i+1, len(splits)):
            rel = arglib.split_relation(splits[i], splits[j])
            if rel == "disjoint":
                compat[-1].append(1)
            elif rel in ("parent", "child", "equal"):
                compat[-1].append(1)
            elif rel == "conflict":
                compat[-1].append(2)
            else:
                raise Exception("unknown rel" + rel)

    # remove redudant conflicts
    for i in xrange(len(splits)):
        for j in xrange(i+1, len(splits)):
            followl = 0
            followr = 0
            dist = j - i
            for k in xrange(1, dist):
                if compat[i][i+k] == 2:
                    followl = 1 # follow left
                if compat[i+k][j] == 2:
                    followr = 1 # follow right
                if followl + followr == 2:
                    compat[i][j] = 0
                    break


    import summon.matrix
    v = summon.matrix.show_dmat(compat, rpart=part, cpart=part, style="quads")


# visualize conflict graph
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 100000      # length of locus
    k = 60         # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen
    mu = 2.5e-8    # mut/site/gen
    u = mu * l     # mut/locus/gen

    arg = arglib.sample_arg(k, n, rho, 0, l)
    mut = arglib.sample_mutations(arg, u)
    mut.sort(key=lambda x: x[2])
    print "sim", len(mut)

    mat = []
    last = None
    lookup = {}
    for tree in arglib.iter_marginal_trees(arg):
        splits = set(arglib.iter_tree_splits(tree))
        for split in splits:
            if split not in lookup:
                lookup[split] = len(lookup)

        if last:
            remove = last - splits
            add = splits - last
            for r in remove:
                for a in add:
                    mat.append((lookup[r], lookup[a], 1))

        last = set(splits)

    v = summon.matrix.show_imat(mat)


# visualize supporting mutations
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 100000      # length of locus
    k = 60         # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen
    mu = 2.5e-8    # mut/site/gen
    u = mu * l     # mut/locus/gen

    arg = arglib.sample_arg(k, n, rho, 0, l)
    mut = arglib.sample_mutations(arg, u)
    mut.sort(key=lambda x: x[2])
    print "sim", len(mut)

    # remove uninformative muts (singletons)
    splits, mut = zip(* (a for a in ((arglib.get_mutation_split(arg, m), m)
                                  for m in mut)
                         if len(a[0]) > 1))


    split_loc = defaultdict(lambda:[])
    for m in mut:
        split = arglib.get_mutation_split(arg, m)
        split_loc[split].append(m[2])

    # get recomb points
    rpos = arglib.get_recomb_pos(arg)

    import summon
    from summon.core import *
    win = summon.Window()

    for i, r in enumerate(rpos):
        splits1 = set(arglib.iter_tree_splits(arg.get_marginal_tree(r-.01)))
        splits2 = set(arglib.iter_tree_splits(arg.get_marginal_tree(r+.01)))
        remove = splits1 - splits2
        add = splits2 - splits1
        all = list(add | remove)

        x = [loc for split in all for loc in split_loc[split]]

        if len(x) > 0:
            win.add_group(lines(color(1,1,1),
                                min(min(x), r), i, max(max(x), r), i))
            for mx in x:
                win.add_group(arglib.draw_mark(mx, i, col=(0,0,1)))

        win.add_group(arglib.draw_mark(r, i, col=(1,0,0)))

    win.home("exact")
