import sys
import util
import copy

from rasmus.progress import Progress

class Smat:
   nrows = 0
   ncols = 0
   nnz = 0
   rows = []
   cols = []
   vals = []
   
   def __init__(self, filename = ""):
      rows = []
      cols = []
      vals = []   
      if (filename != ""):
         self.read(filename)
      
   
   def old__repr__(self):
      string = str(self.nrows) + " " + \
               str(self.ncols) + " " + \
               str(self.nnz) + "\n"
      for i in range(len(self.rows)):         
         string += str(self.rows[i]) + " " + \
                   str(self.cols[i]) + " " + \
                   str(self.vals[i]) + "\n"
      return string

   def read(self, filename):
      infile = file(filename)
      # read header
      (self.nrows, self.ncols, self.nnz) = map(int, infile.next().split())
      
      # read non-zeros
      for line in infile:
         (row, col, val) = line.split()
         self.rows.append(int(row))
         self.cols.append(int(col))
         self.vals.append(float(val))

class AdjMatrix:
    nrows = 0
    ncols = 0
    nnz = 0
    rows = []

    def __init__(self, filename = ""):
       rows = []
       if (filename != ""):
          self.read(filename)


    def __repr__(self):
       string = str(self.nrows) + " " + \
                str(self.ncols) + " " + \
                str(self.nnz) + "\n"
       for x in range(len(self.rows)):
          for y in self.rows[x].keys():
             string += ('%d %d %f\n' % x, y, self.rows[x][y])
       return string
    
    def setDim(self, nrows, ncols, nnz):
        self.nrows = nrows
        self.ncols = ncols
        self.nnz = nnz
        
        # allocate rows array
        for i in range(self.nrows):
            self.rows.append({})

    def read(self, filename, stream = sys.stderr, step = .1):
       infile = file(filename)
       # read header
       (nrows, ncols, nnz) = map(int, infile.next().split())

       self.setDim(nrows, ncols, nnz)

       progress = Progress(0, 0, self.nnz, .01)

       # read non-zeros
       for line in infile:
          progress.update()

          (row, col, val) = line.split()
          if int(row) in range(len(self.rows)):
             self.rows[int(row)][int(col)] = float(val)
          else:
             print "error: %d %s" % (len(self.rows), line)
    
    
    def writeCluto(self, filename, square=False):
       out = util.openStream(filename, "w")

       # write header
       if square:
           out.write("%d %d\n" % (self.nrows, self.nnz))
       else:
           out.write("%d %d %d\n" % (self.nrows, self.ncols, self.nnz))
            

       # write non-zeros in compressed row format
       for i in range(len(self.rows)):
          keys = self.rows[i].keys()
          keys.sort()
          for j in keys:
             out.write("%d %f " % (j+1, self.rows[i][j]))
          out.write("\n")

    def get(self, i, j):
       if j in self.rows[i].keys():
          return self.rows[j][j]
       else:
          return 0
    
    
    def set(self, i, j, value):
        self.rows[i][j] = value
            
    

class Matrix:
    def __init__(self, nrows=0, ncols=0, val=0):
        self.data = []
        self.nrows = nrows
        self.ncols = ncols
        for i in range(nrows):
            self.data.append([])
            for j in range(ncols):
                self.data[i].append(val)
    
    def __getitem__(self, i):
        return self.data[i]
    
    def __repr__(self):
        string = "<MATRIX " + str(self.nrows) + "x" + str(self.ncols) + "\n"
        for row in self.data:
            string += row.__repr__() + "\n"
        string += ">"
        return string

def readLabelMatrix(filename, delim=None, default=0):
    infile = util.openStream(filename)
    mat = util.Dict(2, default)
    
    for line in infile:
        tokens = line.rstrip().split(delim)
        mat[tokens[0]][tokens[1]] = float(tokens[2])
    return mat
    
    
def makeSymetric(mat):
    keys = []
    for i in mat:
        for j in mat[i]:
            keys.append((i,j))
    for i,j in keys:
        mat[j][i] = mat[i][j]
    return mat


def makeMatrix(nrows, ncols, val = 0):
    mat = []
    for i in xrange(nrows):
        row = []
        mat.append(row)
        for j in xrange(ncols):
            row.append(copy.copy(val))
    return mat


def transpose(mat):
    """
    Transpose a matrix
    
    Works better than zip() in that rows are lists not tuples
    """
    
    assert equal(* map(len, mat)), "rows are not equal length"
    
    mat2 = []
    
    for j in xrange(len(mat[0])):
        row2 = []
        mat2.append(row2)
        for row in mat:
            row2.append(row[j])
    
    return mat2


def submatrix(mat, rows=None, cols=None):
    if rows == None:
        rows = range(len(mat))
    if cols == None:
        cols = range(len(mat[0]))
    
    mat2 = []
    
    for i in rows:
        mat2.append([])
        for j in cols:
            mat2[-1].append(mat[i][j])
    
    return mat2


