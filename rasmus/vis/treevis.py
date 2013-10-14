"""
    
    Visualization for phylogenetic trees

"""

# rasmus libs
from rasmus import util
from rasmus import treelib
from compbio import phylo

# summon libs
from summon.core import *
import summon
from summon import sumtree
from summon import shapes
from summon import hud


class TreeViewer (sumtree.SumTree):
    def __init__(self, tree, stree=None, gene2species=None, recon=None,
                 dup_color=(1, 0, 0), loss_color=(0, 0, 1), 
                 **options):
        sumtree.SumTree.__init__(self, tree, **options)
        
        self.mode = "gene"
        self.stree = stree                  # species tree
        self.gene2species = gene2species    # gene to species mapping
        self.recon = recon
        self.bar = None
        
        # colors
        self.dup_color = dup_color
        self.loss_color = loss_color
        
        self.setup_recon(recon)
    
    
    def set_tree(self, tree):
        sumtree.SumTree.set_tree(self, tree)
        self.setup_recon()
                    
    
    def setup_recon(self, recon=None):
        # construct default reconciliation
        if recon == None and self.stree and self.gene2species:
            self.recon = phylo.reconcile(self.tree, self.stree, self.gene2species)
        else:
            self.recon = recon
        
        # construct events
        if self.recon:
            self.events = phylo.label_events(self.tree, self.recon)
            self.losses = phylo.find_loss(self.tree, self.stree, self.recon)
        else:
            self.events = None
            self.losses = None
    
        
    def show(self):
        sumtree.SumTree.show(self)
        
        # set bindings
        self.win.set_binding(input_key("g"), lambda: self.set_mode("gene"))
        self.win.set_binding(input_key("e"), lambda: self.set_mode("events"))
        self.win.set_binding(input_key("r"), lambda: self.set_mode("reroot"))
        self.win.set_binding(input_key("s"), lambda: self.set_mode("swap"))
        self.win.set_binding(input_key("S", "shift"), lambda: self.swap(self.tree.root))
        
        # build sidebar menu
        if self.bar is None:
            self.bar = hud.SideBar(self.win, width=150)
            self.bar.add_item(hud.MenuItem("gene mode (g)", lambda: self.set_mode("gene")))
            self.bar.add_item(hud.MenuItem("events mode (e)", lambda: self.set_mode("events")))
            self.bar.add_item(hud.MenuItem("reroot mode (r)", lambda: self.set_mode("reroot")))
            self.bar.add_item(hud.MenuItem("swap mode (s,S)", lambda: self.set_mode("swap")))
        
        if self.events:
            self.win.add_group(self.draw_events())
    
    
    def draw_events(self):

        # draw duplications
        dups = [color(*self.dup_color)]        
        for node in self.tree:
            if self.events[node] == "dup":
                dups.append(
                    zoom_clamp(
                        shapes.box(node.x - .5, node.y - .5,
                                   node.x + .5, node.y + .5),
                        link=True, link_type="smaller",
                        maxx=8, minx=1,
                        maxy=8, miny=1,
                        origin=(node.x, node.y),
                        prezoom=(self.xscale, 1.0)))
        
        # draw losses
        losses_per_branch = util.hist_dict([node for node, schild in self.losses])
        
        losses = [color(*self.loss_color)]
        for node, nlosses in losses_per_branch.iteritems():
            if node.parent == None:
                continue
                
            x1 = node.parent.x        
            x2 = node.x
            step = (x2 - x1) / float(nlosses + 1)
            
            for x in util.frange(x1 + step, x2-(step/2.0), step):
                losses.append(lines(x, node.y - .2, x, node.y + .2))
        
        return group(group(*dups), group(*losses))


    def set_mode(self, mode):
        print "mode is now '%s'" % mode
        self.mode = mode
        

    def node_click(self, node):
        """Node click callback"""
        
        if self.mode == "gene":
            # print gene names
            #sumtree.SumTree.node_click(self, node)
            print "node: %s\t%f" % (str(node.name), node.dist)
            for key, val in node.data.iteritems():
                print "%s:\t%s" % (key, str(val))
            print
            
        elif self.mode == "events":
            self.print_events(node)
        elif self.mode == "reroot":
            self.reroot(node)
        elif self.mode == "swap":
            self.swap(node)
    
    
    def swap(self, node):
        if len(node.children) > 1:
            node.children = node.children[1:] + [node.children[0]]
        self.set_tree(self.tree)
        self.show()
        self.on_reorder_leaves()        
        
    
    def print_events(self, node):
        """Prints the events that occur on a node"""
        
        print "-" * 20
        print "%s (%.3f sub/site):" % (node.name, node.dist),
        
        if self.events:
            if self.events[node] == "dup":
                print "duplication"
            elif self.events[node] == "spec":
                print "speciation"
            elif self.events[node] == "gene":
                print "extant"
            
            losses = util.groupby(lambda x: x[0], self.losses)
            
            if node in losses:
                print len(losses[node]), "gene losses"
                for i, (node, schild) in enumerate(losses[node]):
                    print "loss %d in species %s" % (i+1, str(schild.name))
                    
                    for sleaf in schild.leaf_names():
                        print "  %s" % sleaf
            
        else:
            print "no reconciliation available"
            
        
    def reroot(self, node):
        try:
            treelib.reroot(self.tree, node.name, newCopy=False)
        except e:
            print e
        print "rerooted on node", node.name
        self.set_tree(self.tree)
        
        self.show()
        self.on_reorder_leaves()
    
    
    def on_reorder_leaves(self):
        """callback for when leaves are reordered"""
        pass
    
    

def read_tree_color_map(filename):
    infile = util.open_stream(filename)
    maps = []
    
    for line in infile:
        expr, red, green, blue = line.rstrip().split("\t")
        maps.append([expr, map(float, (red, green, blue))])
    
    name2color = phylo.make_gene2species(maps)
    
    def leafmap(node):
        return name2color(node.name)

    return treelib.tree_color_map(leafmap)    
readTreeColorMap = read_tree_color_map
