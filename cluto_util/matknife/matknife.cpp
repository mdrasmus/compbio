#include <assert.h>
#include <fstream>
#include <iostream>
#include <list>
#include <vector>
#include <time.h>
#include "HashTable.h"



enum { GRAPH_SIZE = 100000 };
enum { OUT_DEGREE = 10 };
enum {
    MODE_NONE,
    MODE_TO_ADJ,
};

enum {
    MATRIX_NONE,
    MATRIX_ADJ,
    MATRIX_IMAT
};


const char *USAGE = "\
usage: matknife [options]\n\
  -labeled \n\
  -square \n\
  -adj <matrix file>\n\
  -imat <matrix file>\n\
  -to-adj <output matrix file>\n\
  -out-label <out labels file>\n\
  -minval <min value>\n\
\n\
";




class Param
{
public:
    Param() :
        labeled(false),
        square(false),
        minvalue(-1e1000),
        mode(MODE_NONE),
        matType(MATRIX_NONE)
    {}

    bool Parse(int argc, char **argv)
    {
        for (int i=1; i<argc; i++) {
            if (!strcmp(argv[i], "-labeled")) {
                labeled = true;
            } else 
            if (!strcmp(argv[i], "-square")) {
                square = true;
            } else 
            if (!strcmp(argv[i], "-adj")) {
                matrixFile = argv[++i];
                matType = MATRIX_ADJ;
            } else 
            if (!strcmp(argv[i], "-imat")) {
                matrixFile = argv[++i];
                matType = MATRIX_IMAT;
            } else 
            if (!strcmp(argv[i], "-to-adj")) {
                mode = MODE_TO_ADJ;
                outfile = argv[++i];
            } else
            if (!strcmp(argv[i], "-out-label")) {
                outlabelfile = argv[++i];
            } else
            if (!strcmp(argv[i], "-minval")) {
                minvalue = atof(argv[++i]);
            } else {
                fprintf(stderr, "error: unknown option %s\n", argv[i]);
                return false;
            }
        }
        return true;
    }
    
    bool labeled;
    bool square;
    float minvalue;
    int mode;
    int matType;
    string matrixFile;
    string labelFile;
    string outfile;
    string outlabelfile;
};


void WriteStrings(string filename, vector<string> &strs)
{
    FILE *out = fopen(filename.c_str(), "w");
    for (unsigned int i=0; i<strs.size(); i++) {
        fprintf(out, "%s\n", strs[i].c_str());
    }
    fclose(out);
}

bool ReadStrings(string filename, vector<string> &labels)
{
    // open file
    FILE *in;        
    if ((in = fopen(filename.c_str(), "r")) == NULL) {
        fprintf(stderr, "error: cannot open file '%s'\n", filename.c_str());
        return false;
    }
    
    enum { MAXLINE = 1000 };
    char *line = new char [MAXLINE];
    
    while (fgets(line, MAXLINE, in)) {
        int len = strlen(line);
        if (line[len-1] == '\n') 
            line[len-1] = '\0';   // chomp
        labels.push_back(string(line));
    }
    
    delete [] line;
    
    return true;
}

bool ReadInts(string filename, vector<int> &ints)
{
    // open file
    FILE *in;        
    if ((in = fopen(filename.c_str(), "r")) == NULL) {
        fprintf(stderr, "error: cannot open file '%s'\n", filename.c_str());
        return false;
    }
    
    int i;
    while (fscanf(in, "%d", &i) == 1)
        ints.push_back(i);
    
    return true;
}


typedef HashTable<int, float, HashInt> IndexRow;
typedef HashTable<string, int, HashString> Lookup;

class IndexMatrix
{
public:
    IndexMatrix() :
        square(false),
        nedges(0),
        rlookup(GRAPH_SIZE, -1),
        clookup(GRAPH_SIZE, -1)
    {}
    
    inline int AddLabel(Lookup &lookup, vector<string> &labels, string label)
    {
        int id = labels.size();
        lookup[label] = id;
        labels.push_back(label);
        matrix.push_back(new IndexRow);
        return id;
    }
    
    inline int GetId(Lookup &lookup, vector<string> &labels, string label) {
        int id = lookup[label];
        if (id == -1)
            id = AddLabel(lookup, labels, label);
        return id;
    }

    
    bool Read(Param &param, string filename)
    {
        ifstream infile(filename.c_str());
        string label1, label2;
        int id1, id2;
        float value;
        square = param.square;

        printf("reading %s...\n", filename.c_str());

        while (infile.good()) {
            if (param.labeled) {
                infile >> label1 >> label2 >> value;
                id1 = GetId(rlookup, rlabels, label1);
                if (square)
                    id2 = GetId(rlookup, rlabels, label2);
                else
                    id2 = GetId(clookup, clabels, label2);
            } else {
                infile >> label1 >> label2 >> value;
                if (sscanf(label1.c_str(), "%d", &id1) != 1) {
                    fprintf(stderr, "'%s' is not a valid matrix row id", 
                            label1.c_str());
                    return false;
                }
                if (sscanf(label2.c_str(), "%d", &id2) != 1) {
                    fprintf(stderr, "'%s' is not a valid matrix col id", 
                            label2.c_str());
                    return false;
                }

            }
            

            if (value >= param.minvalue) {
                matrix[id1]->Insert(id2, value);
            }
        }

        infile.close();

        return true;
    }
    
    
    void Write(FILE *out)
    {
        for (int i=0; i<NumRows(); i++) {
            IndexRow *row = matrix[i];
            
            for (IndexRow::Iterator j=row->Begin(); j.HasMore();) {
                IndexRow::Node node = j.Next();
                if (square)
                    printf("%s %s %f\n", 
                       rlabels[i].c_str(), 
                       rlabels[node.key].c_str(),
                       node.value);
                else
                    printf("%s %s %f\n", 
                           rlabels[i].c_str(), 
                           clabels[node.key].c_str(),
                           node.value);
            }
        }
    }
    
    
    inline int NumRows()
    { return rlabels.size(); }
    
    inline int NumCols()
    {
        if (square)
            return rlabels.size();
        else
            return clabels.size();
    }
    
    int NumNonZeros()
    {
        int nnz = 0;
        for (unsigned int i=0; i<matrix.size(); i++)
            nnz += matrix[i]->Size();
        return nnz;
    }
    
    bool square;
    int nedges;
    vector<string> rlabels;
    vector<string> clabels;
    Lookup rlookup;
    Lookup clookup;
    vector<IndexRow*> matrix;
};




class AdjMatrix
{
public:
    AdjMatrix() :
        xadj(NULL),
        adjncy(NULL),
        adjwgt(NULL)
    {}
    
    ~AdjMatrix()
    {
        if (xadj)
            delete [] xadj;
        if (adjncy)
            delete [] adjncy;
        if (adjwgt)
            delete [] adjwgt;
    }
    
    
    bool Read(string filename)
    {
        // open file
        FILE *in;        
        if ((in = fopen(filename.c_str(), "r")) == NULL) {
            fprintf(stderr, "error: cannot open file '%s'\n", filename.c_str());
            return false;
        }
        
        enum { MAXLINE = 10000000 };
        char *line = new char [MAXLINE];
        char *delim = " \t";

        
        // read header
        fgets(line, MAXLINE, in);
        int a, b, c;
        if (sscanf(line, "%d %d %d", &a, &b, &c) == 3) {
            nrows = a;
            ncols = b;
            nnz = c;
        } else if (sscanf(line, "%d %d", &a, &b) == 2) {
            nrows = ncols = a;
            nnz = b;
        } else {
            fprintf(stderr, "bad format on line 1 '%s'\n", filename.c_str());
            delete [] line;
            return false;
        }
        
        printf("reading matrix '%s' nrows=%d, ncols=%d, nnz=%d...\n", 
               filename.c_str(), nrows, ncols, nnz);
        
        // allocate space
        xadj = new int [nrows + 1];
        adjncy = new int [nnz];
        adjwgt = new float [nnz];
        
        // read data
        int i = 0;
        int j = 0;
        while (fgets(line, MAXLINE, in)) {
            char *ptr = line;
            char *token1, *token2;
            xadj[i++] = j;
            
            // skip lead white space
            while (ptr[0] == ' ' || ptr[0] == '\t')
                strsep(&ptr, delim);
            
            while ((token1 = strsep(&ptr, delim)) &&
                   (token2 = strsep(&ptr, delim)))
            {
                int col = atoi(token1) - 1;
                assert(col >= 0 && col < ncols);
                adjncy[j] = col;
                adjwgt[j] = atof(token2);
                j++;
            }
        }
        assert(j == nnz);
        assert(i == nrows || i == nrows+1);
        xadj[nrows] = nnz;
        
        // clean up
        delete [] line;
        fclose(in);
        
        return true;
    }
    
    
    void Submatrix(AdjMatrix &matrix, vector<int> &verts)
    {
        // create bool array for verts
        int *selected = new int [matrix.nrows];
        for (int i=0; i<matrix.nrows; i++)
            selected[i] = -1;
        for (unsigned int i=0; i<verts.size(); i++)
            selected[verts[i]] = i;
        
        // count new nnz
        nrows = verts.size();        
        nnz = 0;
        for (int i=0; i<nrows; i++)
            for (int j=matrix.xadj[verts[i]]; j<matrix.xadj[verts[i]+1]; j++)
                if (selected[matrix.adjncy[j]] != -1)
                    nnz++;
        
        // allocate space
        xadj = new int [nrows + 1];
        adjncy = new int [nnz];
        adjwgt = new float [nnz];
        
        // copy over the submatrix
        int k=0;
        for (int i=0; i<nrows; i++) {
            xadj[i] = k;
            for (int j=matrix.xadj[verts[i]]; j<matrix.xadj[verts[i]+1]; j++) {
                if (selected[matrix.adjncy[j]] != -1) {
                    adjncy[k] = selected[matrix.adjncy[j]];
                    adjwgt[k] = matrix.adjwgt[j];
                    k++;
                }
            }
        }
        assert(k == nnz);
        xadj[nrows] = nnz;
        
        // clean up
        delete [] selected;
    }
    
    void Write(FILE *out)
    {
        fprintf(out, "%d %d", nrows, nnz);
        for (int i=0; i<nrows; i++) {
            for (int j=xadj[i]; j<xadj[i+1]; j++) {
                fprintf(out, "%d %f ", adjncy[j]+1, adjwgt[j]);
            }
            fprintf(out, "\n");
        }
    }
    
    int nrows;
    int ncols;
    int nnz;    
    int *xadj;
    int *adjncy;
    float *adjwgt;
};


void Convert(Param &param, IndexMatrix &imatrix, string outfile, string labelfile)
{
    // count nnz
    int nnz = imatrix.NumNonZeros();
    
    FILE *out = fopen(outfile.c_str(), "w");
    
    // write matrix
    if (imatrix.square)
        fprintf(out, "%d %d\n", imatrix.NumRows(), nnz);    
    else
        fprintf(out, "%d %d %d\n", imatrix.NumRows(), imatrix.NumCols(), nnz);    
    for (int i=0; i<imatrix.NumRows(); i++) {        
        for (IndexRow::Iterator j=imatrix.matrix[i]->Begin(); j.HasMore();) {
            IndexRow::Node node = j.Next();
            fprintf(out, "%d %f ", node.key + 1, node.value);
        }
        fprintf(out, "\n");
    }
    
    fclose(out);
    
    // write labels
    if (labelfile != "") {

        if (imatrix.square) {
            WriteStrings(labelfile, imatrix.rlabels);
        } else {
            string rlabelfile = labelfile + ".rlabel";
            string clabelfile = labelfile + ".clabel";                        
            
            WriteStrings(rlabelfile, imatrix.rlabels);
            WriteStrings(clabelfile, imatrix.clabels);
        }
        
        
    }
}







int main(int argc, char **argv)
{
    // read parameters
    Param param;
    if (!param.Parse(argc, argv)) {
        fprintf(stderr, USAGE);
        return 1;
    }
    
    if (param.mode == MODE_NONE) {
        fprintf(stderr, USAGE);
        return 1;
    }
    
    // load matrix
    IndexMatrix imatrix;
    AdjMatrix amatrix;
    
    switch (param.matType) {
        case MATRIX_ADJ:
            if (!amatrix.Read(param.matrixFile))
                exit(1);
            break;
        case MATRIX_IMAT:
            if (!imatrix.Read(param, param.matrixFile))
                exit(1);
            break;
        default:
            fprintf(stderr, USAGE);
            return 1;
    }
    
    // execute mode
    switch (param.mode) {
        case MODE_TO_ADJ:
            if (param.matType != MATRIX_IMAT) {
                fprintf(stderr, "must input index matrix (-imat)\n");
                exit(1);
            }
            Convert(param, imatrix, param.outfile, param.outlabelfile);
        break;
    }
}

