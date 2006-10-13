import util

def repeatedDfs(vertices, visit, getNeighbors, newComponent = lambda x: 0):
    openset = []
    closedset = {}
    component = -1

    for start in vertices:
        # skip closed vertices
        if start in closedset:
            continue
        
        component += 1
        newComponent(component)
        openset = [start]
        dfs(vertices, visit, getNeighbors, openset, closedset)
        

def dfs(vertices, visit, getNeighbors, openset, closedset = None):
    """depth first search"""
    
    if closedset == None:
        closedset = {}
    
    # search from openset
    while len(openset) > 0:
        vertex = openset.pop()

        # skip closed vertices
        if vertex in closedset:
            continue

        # visit vertex
        visit(vertex)

        openset.extend(getNeighbors(vertex))

        # close new vertex
        closedset[vertex] = 1


def bfs(vertices, visit, getNeighbors, openset, closedset = None):
    """breadth first search"""
    
    if closedset == None:
        closedset = {}
    
    # search from openset
    while len(openset) > 0:
        print "openset", openset
        vertex = openset[0]
        openset = openset[1:]

        # skip closed vertices
        if vertex in closedset:
            continue

        # visit vertex
        visit(vertex)

        openset.extend(getNeighbors(vertex))
        
        # close new vertex
        closedset[vertex] = 1
        

def connectedComponents(vertices, neighborFunc):
    components = []
    
    def newComponent(component):
        components.append([])
    
    def visit(vertex):
        components[-1].append(vertex)

    repeatedDfs(vertices, visit, neighborFunc, newComponent)
    
    return components


def neighborhood(vertex, getNeighbors, radius):
    openset = [(vertex , 0)]
    closedset = {}
    
    while len(openset) > 0:
        (v, dist) = openset.pop()
        
        # stop if v is already visited
        if v in closedset and closedset[v] <= dist:
            continue
        
        if dist < radius:
            next = getNeighbors(v)
            for i in next:
                openset.append((i, dist+1))
        
        closedset[v] = dist
    return closedset


def subgraph(vertices, outEdges):
    # create empty matrix
    mat = {}
    vertices = set(vertices)
        
    for v in vertices:
        mat[v] = {}
        edges = outEdges(v)
        for (u, edge) in edges:
            if u in vertices:
                mat[v][u] = edge
    
    return mat


if __name__ == "__main__":
    g = {1: [2, 3], 
         2: [1, 3],
         3: [1, 2],
         4: [5, 6],
         5: [4, 6],
         6: [4, 5],
         7: [8, 9],
         8: [7, 10, 11],
         9: [7, 12, 13],
         10: [8],
         11: [8],
         12: [9],
         13: [9]}
    
    print connectedComponents(g, lambda x: g[x])
    
    def visit(x):
        print "visit", x
    
    bfs(g, visit, lambda x: g[x], [7])

    

