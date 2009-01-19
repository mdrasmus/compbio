# python libs
import math
import random
import sys



#=============================================================================

class UnionFind:
    """An implementation of the UNINON/FIND algorithm"""

    def __init__(self, items):
        self.parent = None    
        self.items = dict.fromkeys(items, 1)

    def __contains__(self):
        return item in self.root().items

    def __len__(self):
        return len(self.root().items)
    
    def __iter__(self):
        return iter(self.root().items)
    
    
    def add(self, item):
        self.root().items[item] = 1
    
    def root(self):
        node = self
        while node.parent:
            node = node.parent
        if node != self:
            self.parent = node
        return node
    
    def same(self, other):
        return self.root() == other.root()
    
    def union(self, other):
        root1 = self.root()
        root2 = other.root()
        if root1 == root2:
            return
        
        root1.items.update(root2.items)
        root2.items = {}
        root2.parent = root1
    
    def members(self):
        return self.root().items.keys()
    
    
    # old function DON'T USE
    
    def has(self, item):
        """DEPRECATED: use x in set"""
        return item in self.members()
    
    def size(self):
        """DEPRECATED: use len(set)"""
        return len(self.root().items)


#=============================================================================
# QuadTree data structure
    
class Rect:
    """A representation of a rectangle"""       
    
    def __init__(self, x1, y1, x2, y2):
        if x1 < x2:
            self.x1 = x1
            self.x2 = x2
        else:
            self.x1 = x2
            self.x2 = x1
        if y1 < y2:
            self.y1 = y1
            self.y2 = y2
        else:
            self.y1 = y2
            self.y2 = y1

class QuadNode:
    item = None
    rect = None
    
    def __init__(self, item, rect):
        self.item = item
        self.rect = rect
        
        
class QuadTree:
    MAX = 10
    MAX_DEPTH = 10
    
    def __init__(self, x, y, size, depth = 0):
        self.nodes = []
        self.children = []
        self.center = [x, y]
        self.size = size
        self.depth = depth
    
    def insert(self, item, rect):
        if len(self.children) == 0:
            self.nodes.append(QuadNode(item, rect))
            
            if len(self.nodes) > self.MAX and self.depth < self.MAX_DEPTH:
                self.split()
        else:
            self.insertIntoChildren(item, rect)
    
    def insertIntoChildren(self, item, rect):
        if rect.x1 < self.center[0]:
            if rect.y1 < self.center[1]:
                self.children[0].insert(item, rect)
            if rect.y2 > self.center[1]:
                self.children[1].insert(item, rect)
        if rect.x2 > self.center[0]:
            if rect.y1 < self.center[1]:
                self.children[2].insert(item, rect)
            if rect.y2 > self.center[1]:
                self.children[3].insert(item, rect)
                   
    def split(self):
        self.children = [QuadTree(self.center[0] - self.size/2,
                                  self.center[1] - self.size/2,
                                  self.size/2, self.depth + 1),
                         QuadTree(self.center[0] - self.size/2,
                                  self.center[1] + self.size/2,
                                  self.size/2, self.depth + 1),
                         QuadTree(self.center[0] + self.size/2,
                                  self.center[1] - self.size/2,
                                  self.size/2, self.depth + 1),
                         QuadTree(self.center[0] + self.size/2,
                                  self.center[1] + self.size/2,
                                  self.size/2, self.depth + 1)]

        for node in self.nodes:
            self.insertIntoChildren(node.item, node.rect)
        self.nodes = []

    def query(self, rect, results = {}, ret = True):
        if ret:
            results = {}
        
        if len(self.children) > 0:
            if rect.x1 < self.center[0]:
                if rect.y1 < self.center[1]:
                    self.children[0].query(rect, results, False)
                if rect.y2 > self.center[1]:
                    self.children[1].query(rect, results, False)
            if rect.x2 > self.center[0]:
                if rect.y1 < self.center[1]:
                    self.children[2].query(rect, results, False)
                if rect.y2 > self.center[1]:
                    self.children[3].query(rect, results, False)
        else:
            for node in self.nodes:
                if node.rect.x2 > rect.x1 and node.rect.x1 < rect.x2 and \
                   node.rect.y2 > rect.y1 and node.rect.y1 < rect.y2:
                    results[node.item] = True
                    
        if ret:
            return results.keys()
                    
    def getSize(self):
        size = 0
        for child in self.children:
            size += child.getSize()
        size += len(self.nodes)
        return size





#=============================================================================
# TODO: make a funtion based linear search

def binsearch(lst, val, compare=cmp, order=1):
    """Performs binary search for val in lst using compare
    
       if val in lst:
          Returns (i, i) where lst[i] == val
       if val not in lst  
          Returns index i,j where
            lst[i] < val < lst[j]
        
       runs in O(log n)
    """
    
    assert order == 1 or order == -1
    
    low = 0
    top = len(lst) - 1
    
    if len(lst) == 0:
        return None, None
    
    if compare(lst[-1], val) * order == -1:
        return (top, None)
    
    if compare(lst[0], val) * order == 1:
        return (None, low)
    
    while top - low > 1:
        ptr = (top + low) / 2
        
        comp = compare(lst[ptr], val) * order
        
        if comp == 0:
            # have we found val exactly?
            return ptr, ptr
        elif comp == -1:
            # is val above ptr?
            low = ptr
        else:
            top = ptr
            
    
    # check top and low for exact hits
    if compare(lst[low], val) == 0:
        return low, low
    elif compare(lst[top], val) == 0:
        return top, top
    else:
        return low, top



if __name__ == "__main__":
    
    if True:
        set1 = UnionFind()
        set2 = UnionFind()
        set3 = UnionFind()

        set1.add(1)
        set1.add(2)
        print set1.size()
        set2.add(3)
        set2.add(4)
        set2.add(5)    
        print set2.size()
        set3.add(5)
        set3.add(6)
        set3.add(7)
        print set3.size()    
        print set1.same(set2)    
        set1.union(set2)
        print set1.same(set2)
        set1.union(set3)

        print set1.members()
        print set1.size(), set2.size()

