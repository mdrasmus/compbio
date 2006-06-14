/*******************************************************************************
 SINDER - Matt Rasmussen
 Tree.h
 11/19/05

*******************************************************************************/

#include <vector>
#include "Tree.h"
#include "Matrix.h"
#include "common.h"


namespace Sinder {

using namespace std;


float TreeNode::Height()
{
    float height = 0.0;
    
    for (unsigned int i=0; i<m_children.size(); i++) {
        float height2 = m_children[i]->Height();
        if (height2 > height)
            height = height2;
    }
    
    return height + m_dist;
}



TreeNode *Tree::AddNode(TreeNode* node=NULL, TreeNode *parent=NULL)
{
    if (!node)
        node = new TreeNode();
    
    node->SetId(m_nodes.size());
    m_nodes.push_back(node);
    
    if (parent) {
        parent->AddChild(node);
    }
    return node;
}


int GetSizesHelper(TreeNode *node, int *sizes)
{
    int size = 0;
    
    for (int i=0; i<node->GetNumChildren(); i++)
        size += GetSizesHelper(node->GetChild(i), sizes);
    
    // include self in size and return
    size++;
    sizes[node->GetId()] = size;
    return size;
}


int *Tree::GetSizes()
{
    int *sizes = new int [m_nodes.size()];
    
    GetSizesHelper(m_root, sizes);
    
    return sizes;
}



char ReadChar(FILE *stream, int &depth)
{
    char chr;
    do {
        if (fread(&chr, sizeof(char), 1, stream) != 1) {
            // indicate EOF
            return '\0';
        }
    } while (chr == ' ' && chr == '\n');
    
    // keep track of paren depth
    if (chr == '(') depth++;
    if (chr == ')') depth--;
    
    return chr;
}


char ReadUntil(FILE *stream, string &token, char *stops, int &depth)
{
    char chr;
    token = "";
    while (true) {
        chr = ReadChar(stream, depth);
        if (!chr)
            return chr;
        
        // compare char to stop characters
        for (char *i=stops; *i; i++) {
            if (chr == *i)
                return chr;
        }
        token += chr;
    }
}

void Indent(FILE *stream, int indent)
{
    for (int i=0; i<indent; i++)
        fprintf(stream, " ");
}


string Trim(const char *word)
{
    char buf[101];
    sscanf(word, "%100s", buf);
    return buf;
}


float ReadDist(FILE *infile, int &depth)
{
    float dist = 0;
    fscanf(infile, "%f", &dist);
    return dist;
}


TreeNode *Tree::ReadNode(FILE *infile, TreeNode *parent, int &depth)
{
    char chr, char1;
    TreeNode *node;
    string token;

    // read first character
    if (!(char1  = ReadChar(infile, depth))) {
        Error("unexpected end of file");
        return NULL;
    }
    

    if (char1 == '(') {
        // read internal node
    
        int depth2 = depth;
        node = AddNode(NULL, parent);
        
        // read all child nodes at this depth
        while (depth == depth2) {
            TreeNode *child = ReadNode(infile, node, depth);
            if (!child)
                return NULL;
        }
        
        // read distance for this node
        char chr = ReadUntil(infile, token, "):,", depth);
        if (chr == ':')
            node->SetDistance(ReadDist(infile, depth));
        if (!(chr = ReadUntil(infile, token, "):,", depth)))
            return NULL;
        return node;
    } else {
        // read leaf
        
        node = AddNode(NULL, parent);
        
        // read name
        if (!(chr = ReadUntil(infile, token, ":),", depth)))
            return NULL;
        token = char1 + Trim(token.c_str());
        node->SetName(token);
        
        // read distance for this node
        if (chr == ':')
            node->SetDistance(ReadDist(infile, depth));
        if (!(chr = ReadUntil(infile, token, ":),", depth)))
            return NULL;
        return node;
    }
}



bool Tree::Read(FILE *infile)
{
    TreeNode *node;
    int depth = 0;
    string token;
    
    // init tree with root
    m_root = AddNode();
    
    // ensure that tree begins with open paren
    char chr = ReadUntil(infile, token, "(", depth);    
    if (chr != '(')
        return false;
    
    // add nodes to root
    while ((depth > 0) && (node = ReadNode(infile, m_root, depth)));
    
    // return success status
    return depth == 0;
}


void Tree::WriteNode(FILE *outfile, TreeNode *node, int indent)
{
    // write node
    if (node->IsLeaf()) {
        Indent(outfile, indent);
        fprintf(outfile, "%s", node->GetName().c_str());
    } else {
        // print open paren
        Indent(outfile, indent);
        fprintf(outfile, "(\n");
        
        // write children
        int nchild = node->GetNumChildren();
        for (int i=0; i<nchild; i++) {
            WriteNode(outfile, node->GetChild(i), indent+1);
            if (i < nchild-1)
                fprintf(outfile, ",\n");
        }
        
        // print close paren
        fprintf(outfile, "\n");
        Indent(outfile, indent);
        fprintf(outfile, ")");
    }
    
    // write node's distance
    if (node != m_root) {
        fprintf(outfile, ":%f", node->GetDistance());
    }
}


void Tree::Write(FILE *outfile)
{
    if (m_root) {
        WriteNode(outfile, m_root, 0);
        fprintf(outfile, ";\n");
    }
}


void DrawLine(Matrix<char> &matrix, char chr, int x1, int y1, int x2, int y2)
{
    float stepx, stepy;
    int steps;
    
    steps = max(x2-x1, y2-y1);
    
    stepx = float(x2 - x1+1) / steps;
    stepy = float(y2 - y1+1) / steps;
    
    float x = x1;
    float y = y1;
    for (int i=0; i<steps; i++) {
        matrix[int(y)][int(x)] = chr;
        x += stepx;
        y += stepy;
    }
}



// TODO: finish display
int DisplayNode(Matrix<char> &matrix, TreeNode *node, 
                 int x, int y, float scale)
{
    //int topy = y;
    //int topx = x;
    x += int(scale * node->GetDistance());
    
    for (int i=0; i<node->GetNumChildren(); i++) {
        y += DisplayNode(matrix, node->GetChild(i), x, y, scale);
    }
    
    return 0;
}


// TODO: finish display
void Tree::Display(FILE *outfile, int width, float scale)
{
    Matrix<char> matrix(Size(), width);
    matrix.SetAll(' ');
        
    // set default scaling if needed
    if (scale == -1)
        scale = (width-1) / m_root->Height();
    
    DisplayNode(matrix, m_root, 0, 0, scale);
    
    // write out matrix
    for (int i=0; i<matrix.NumRows(); i++) {
        for (int j=0; j<matrix.NumCols(); j++)
            fprintf(outfile, "%c", matrix[i][j]);
        fprintf(outfile, "\n");
    }
}


};
