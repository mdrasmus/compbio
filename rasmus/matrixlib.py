"""

    Common dense and sparse matrix input/output/conversion

    imat  -- sparse index matrix
             [(i, j, v), ...]
             
    rmat  -- sparse row compressed matrix
             [{j1: v1, j2: v2, ...}, ...]
             
    dmat  -- dense matrix
             [[v11, v12, ...], [v21, v22, ...], ...]
             
    lmat  -- sparse labeled matrix
             {"row1": {"col2": v, ...}, ...}

    dlmat -- dense labeled matrix
             [["", "col1", "col2", ...],
              ["row1", v11, v12, ...],
              ["row2", v21, v22, ...]]
             
    ilmat -- iterate sparse label matrix
             [(labeli, labelj, v), ...]


"""

# python libs
from collections import defaultdict
import copy



def make_matrix(nrows, ncols, val=0):
    """Returns a dense matrix of desired size"""
    
    mat = []
    for i in xrange(nrows):
        row = []
        mat.append(row)
        for j in xrange(ncols):
            row.append(copy.copy(val))
    return mat


def make_lmat(default=0):
    """Returns an empty lmat"""
    return defaultdict(lambda: defaultdict(lambda: default))
    

def submatrix(mat, rows=None, cols=None):
    """Returns the submatrix of mat"""
    
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


def transpose(mat):
    """
    Transpose a matrix
    
    Works better than zip() in that returned rows are lists not tuples
    """
    
    assert len(set(map(len, mat))) == 1, "rows are not equal length"
    
    mat2 = []
    
    for j in xrange(len(mat[0])):
        row2 = []
        mat2.append(row2)
        for row in mat:
            row2.append(row[j])
    
    return mat2


def transpose_imat(nrows, ncols, nnz, imat):
    """Transpose an index matrix iterator"""

    def data():
        for i, j, v in imat:
            yield j, i, v

    return ncols, nrows, nnz, data()




#=============================================================================
# matrix conversion


def imat2rmat(nrows, ncols, nnz, imat):
    """Converts imat iterator to rmat"""

    rows = []
    for i in xrange(nrows):
        rows.append({})

    # store data
    for i, j, v in imat:
        rows[i][j] = v

    return rows

def rmat2imat(nrows, ncols, nnz, rmat):
    """Converts rmat to imat iterator"""

    for i, row in enumerate(rmat):
        for j, v in row.iteritems():
            yield i, j, v


def dmat2imat(dmat):
    """Converts dense matrix to sparse index matrix"""
    
    for i, row in enumerate(dmat):
        for j, v in enumerate(row):
            yield i, j, v

def imat2dmat(nrows, ncols, nnz, imat):
    """Converts sparse index matrix to a dense matrix"""

    dmat = make_matrix(nrows, ncols)
    for i, j, v in imat:
        dmat[i][j] = v

    return dmat


def ilmat2lmat(ilmat, default=0):
    """Converts a labeled matrix iterator (ilmat) to a dict of dicts (lmat)"""
    
    lmat = defaultdict(lambda: defaultdict(lambda: default))
    for r, c, v in ilmat:
        lmat[r][c] = v
    return lmat


def lmat2ilmat(lmat):
    """Converts a dict of dicts (lmat) to labeled matrix iterator (ilmat)"""

    for row in lmat:
        for col, val in lmat[row].iteritems():
            yield row, col, val


def ilmat2imat(ilmat, rowlabels, collabels):
    """
    Converts labeled matrix iterator (ilmat) to indexed matrix iterator (imat)
    """

    rowlookup = dict((l, i) for i, l in enumerate(rowlabels))
    collookup = dict((l, i) for i, l in enumerate(rowlabels))

    for r, c, v in ilmat:
        yield rowlookup[r], collookup[c], v


def lmat2imat(lmat, rowlabels, collabels):
    """
    Converts a dict of dicts (lmat) to a indexed matrix iterator (imat)
    """
    return ilmat2imat(lmat2ilmat(lmat), rowlabels, collabels)


def imat2ilmat(imat, rowlabels, collabels):
    """
    Converts indexed matrix iterator (imat) to labeled matrix iterator (ilmat)
    """

    for i, j, v in imat:
        yield rowlabels[i], collabels[j], v


def imat2dlmat(imat, rowlabels, collabels):
    """
    Converts indexed matrix iterator (imat) to dense labeled matrix (dlmat)
    """

    dlmat = make_matrix(len(rowlabels)+1, len(collabels)+1)
    rowlookup = dict((l, i+1) for i, l in enumerate(rowlabels))
    collookup = dict((l, i+1) for i, l in enumerate(rowlabels))    

    # set labels
    dlmat[0][0] = ""
    for i, label in enumerate(rowlabels):
        dlmat[i+1][0] = label
    for i, label in enumerate(collabels):
        dlmat[0][i+1] = label

    # set values
    for i, j, v in imat:
        dlmat[rowlookup[i]][collookup[j]] = v

    return dlmat


def dlmat2imat(dlmat):
    """
    Converts dense labeled matrix (dlmat) to indexed matrix iterator (imat)
    """

    for i, row in enumerate(dlmat):
        if i > 0:
            for j, v in enumerate(row):
                if j > 0:
                    yield i-1, j-1, v
    

def dmat2dlmat(dmat, rowlabels, collabels):
    """
    Converts dense matrix (dmat) to dense labeled matrix (dlmat)
    """
    
    dlmat = make_matrix(len(rowlabels)+1, len(collabels)+1)
    rowlookup = dict((l, i+1) for i, l in enumerate(rowlabels))
    collookup = dict((l, i+1) for i, l in enumerate(collabels))    

    # set labels
    dlmat[0][0] = ""
    for i, label in enumerate(rowlabels):
        dlmat[i+1][0] = label
    for i, label in enumerate(collabels):
        dlmat[0][i+1] = label

    # set values
    for i in xrange(len(rowlabels)):
        for j in xrange(len(collabels)):
            dlmat[i+1][j+1] = dmat[i][j]

    return dlmat


#=============================================================================
# dense matrix I/O

def parse_dmat_header(first_row, header, rows):

    # parse possible header
    if header is None:
        # auto detect header
        if len(first_row) == 1:
            nrows = int(first_row[0])
            ncols = nrows
        elif len(first_row) == 2:
            if len(first_row) != len(rows[1]):
                nrows = int(first_row[0])
                ncols = int(first_row[1])
            else:
                raise Exception("Cannot automatically detect matrix header")

        else:
            # no header
            # infer matric dimension
            nrows = ncols = None
        
    elif header:
        # explicitly parse header        
        if len(first_row) == 1:
            nrows = int(first_row[0])
            ncols = nrows
        elif len(first_row) == 2:
            nrows = int(first_row[0])
            ncols = int(first_row[1])
        else:
            raise Exception("Wrong number of entries in header (expected 1 or 2)")
            
    else:
        # no header
        # infer matrix dimension
        nrows = ncols = None

    return nrows, ncols


def read_dmat(infile, header=False):
    """
    Reads dense matrix
    header -- can be True,False,None specifying whether to read a header (True)
              or auto-detect it (None)
    """

    # read file
    rows = [line.rstrip().split() for line in infile]    
    nrows, ncols = parse_dmat_header(rows[0], header, rows)
    if nrows is not None:
        # skip header
        rows = rows[1:]
    
    # convert data
    data = [[float(v) for v in row] for row in rows]
    
    # assert that all rows have the same number of values
    assert len(set(map(len, data))) == 1

    # assert dimensions
    if nrows is None:
        nrows = len(data)
    else:
        assert nrows == len(data), "wrong number of rows"

    if ncols is None:
        ncols = len(data[0])
    else:
        assert ncols == len(data[0]), "wrong number of columns"
    
    # return data
    nnz = nrows * ncols
    return nrows, ncols, nnz, data


def iter_dmat(infile, header=False):
    """
    Iterates a dense matrix
    header -- can be True,False,None specifying whether to read a header (True)
              or auto-detect it (None)

    Returns nrows, ncols, nnz, imat (index matrix interator)
    """
    
    # read file
    rows = [line.rstrip().split() for line in infile]    
    nrows, ncols = parse_dmat_header(rows[0], header, rows)
    if nrows is not None:
        # skip header
        nnz = nrows * ncols
        rows = rows[1:]
    else:
        nnz = None
    
    def data():
        for i, row in enumerate(rows):
            for j, v in enumerate(row):
                yield i, j, float(v)

    return nrows, ncols, nnz, data()


def write_dmat(out, dmat, square=False):
    """Write a dense matrix file"""

    # write header
    if square:
        assert len(dmat) == len(dmat[0])
        out.write("%d\n" % len(dmat))
    else:
        out.write("%d\t%d\n" % (len(dmat), len(dmat[0])))

    # write data
    for row in dmat:
        for val in row[:-1]:
            out.write("%f\t" % val)
        out.write("%f\n" % row[-1])


#=============================================================================
# index matrix I/O

def iter_imat(infile):
    """
    Read a sparse index matrix
    Returns nrows, ncols, nnz, imat (index matrix iterator)
    """
    
    line = infile.next()
    tokens = line.rstrip().split()

    try:
        if len(tokens) == 3:
            nrows, ncols, nnz = map(int, tokens)
        elif len(tokens) == 2:
            nrows, nnz = map(int, tokens)
            ncols = nrows
    except:
        raise Exception("header error: expected (nrows, nnz) or (nrows, ncols, nnz) in first line")

    def data():
        for line in infile:
            tokens = line.split()
            yield int(tokens[0]), int(tokens[1]), float(tokens[2])

    return nrows, ncols, nnz, data()


def read_imat(infile):
    """
    Read a sparse index matrix
    Returns nrows, ncols, nnz, (rows, cols, vals)
    """

    nrows, ncols, nnz, imat = iter_imat(infile)
    return nrows, ncols, nnz, list(imat)
    

def write_imat(out, nrows, ncols, nnz, imat):
    """Write an index matrix"""

    out.write("%d\t%d\t%d\n" % (nrows, ncols, nnz))
    for i, j, v in imat:
        out.write("%d\t%d\t%f\n" % (i, j, v))


#=============================================================================
# compressed-row matrix I/O

def iter_rmat(infile):
    """
    Read an compressed-row matrix
    Columns are 1 indexed in file, but 0 index in memory
    Returns nrows, ncols, nnz, imat (index matrix iterator)    
    """

    line = infile.next()
    tokens = map(int, line.rstrip().split())

    if len(tokens) == 2:
        nrows, nnz = tokens
        ncols = nrows
    elif len(tokens) == 3:
        nrows, ncols, nnz = tokens
    else:
        raise Exception("wrong number of fields in header (must 2 or 3)")

    def data():
        for i, line in enumerate(infile):
            tokens = line.split()

            for j in xrange(0, len(tokens)-1, 2):
                yield i, int(tokens[j]) - 1, float(tokens[j+1])

    return nrows, ncols, nnz, data()


def read_rmat(infile):
    """
    Read an compressed-row matrix
    Columns are 1 indexed in file, but 0 index in memory
    Returns nrows, ncols, nnz, rmat (list of dicts)
    """

    nrows, ncols, nnz, data = iter_rmat(infile)
    return nrows, ncols, nnz, imat2rmat(nrows, ncols, nnz, data)


def write_rmat(out, nrows, ncols, nnz, rmat, square=False):
    """
    Write a compressed-row matrix
    Columns are 1 indexed in file, but 0 index in memory
    """

    if square:
        assert nrows == ncols
        out.write("%d\t%d\n" % (nrows, nnz))
    else:
        out.write("%d\t%d\t%d\n" % (nrows, ncols, nnz))

    for row in rmat:
        items = row.items()
        items.sort()

        if len(items) > 0:
            for j, v in items[:-1]:
                out.write("%d\t%f\t" % (j + 1, v))
            out.write("%d\t%f\n" % (items[-1][0] + 1, items[-1][1]))
        else:
            out.write("\n")


#=============================================================================
# labeled sparse matrix I/O

def read_lmat(infile, delim=None, default=0):
    """
    Reads a labeled sparsed matrix
    Returns matrix as a dict of dicts
    """
    return ilmat2lmat(iter_lmat(infile, delim=delim), default=default)


def iter_lmat(infile, delim=None, default=0):
    """
    Read a labeled sparsed matrix
    Returns an labeled matrix iterator, ilmat = (rowlabel, collabel, value)
    """

    for line in infile:
        tokens = line.rstrip().split(delim)
        yield tokens[0], tokens[1], float(tokens[2])


def write_lmat(out, lmat):
    """
    Writes a labeled sparsed matrix
    """
    
    for row in lmat:
        for col, val in lmat[row].iteritems():
            out.write("%s\t%s\t%f\n" % (row, col, val))
            
#=============================================================================
# dense labeled matrix I/O

def read_dlmat(infile, delim=None, default=0):

    infile = iter(infile)

    # read col labels
    dlmat = [infile.next().rstrip().split(delim)]

    for line in infile:
        tokens = line.rstrip().split(delim)
        dlmat = [tokens[0]] + map(float, tokens[1:])

    return dlmat


def write_dlmat(out, dlmat):
    
    # write data
    for row in dlmat:
        for i in xrange(len(row)-1):
            out.write(str(row[i]) + "\t")
        out.write(str(row[-1]) + "\n")

