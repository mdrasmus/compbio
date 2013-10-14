import sys
import os


from rasmus import util



def writeGraphviz(mat, out=sys.stdout, format="undirected", 
                  param="overlap=\"false\";"):
    if format == "undirected":
        print >>out, "graph G {"
        print >>out, param
        for i in mat:
            for j in mat[i]:
                if i < j:
                    print >>out, i, "--", j, ";"
        print >>out, "}"
    elif format == "directed":
        print >>out, "digraph G {"
        for i in mat:
            for j in mat[i]:
               print >>out, i, "->", j, ";"
        print >>out, "}"        

def visualize(mat, outfile, format="undirected", 
              options="-Tjpg", param="overlap=\"false\";"):
    if format == "undirected":
        out = os.popen("neato "+ options + " -o " + outfile, "w")
    elif format == "directed":
        out = os.popen("dot "+ options + " -o " + outfile, "w")
    print out, format
    writeGraphviz(mat, out, format, param)
    

if __name__ == "__main__":
    mat = util.Dict(dim=2)

    mat[1][2] = 1
    mat[1][3] = 1
    mat[2][4] = 1
    mat[3][4] = 1
    mat[4][5] = 1

    #writeGraphviz(mat, sys.stdout, "directed")
    visualize(mat, "out2.jpg", "directed")


class GraphViz:
    pass


