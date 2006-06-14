/*******************************************************************************
 SINDIR - Matt Rasmussen
 Tree.h
 11/19/05

*******************************************************************************/

#include <stdio.h>
#include <string>
#include <vector>




namespace Sinder {

using namespace std;

class TreeNode
{
public:
    TreeNode(int id=0, string name="", float dist=0.0) :
        m_id(id),
        m_name(name),
        m_dist(dist),
        m_parent(NULL)
    {}

    // accessors
    inline int GetId() { return m_id; }
    inline void SetId(int id) { m_id = id; }
    inline string GetName() { return m_name; }
    inline string SetName(string name) { return m_name = name; }
    inline float GetDistance() { return m_dist; }
    inline void SetDistance(float dist) { m_dist = dist; }
    inline int GetNumChildren() { return m_children.size(); }
    inline TreeNode *GetChild(int i) { return m_children[i]; }
    inline void SetChild(int i, TreeNode *child)
    {
        m_children[i] = child; 
        child->SetParent(this);
    }
    inline void AddChild(TreeNode *child)
    {
        m_children.push_back(child);
        child->SetParent(this);
    }
    inline TreeNode *GetParent() { return m_parent; }
    inline void SetParent(TreeNode *parent) { m_parent = parent; }
    inline bool IsLeaf() { return m_children.size() == 0; }
    
    float Height();
    
private:
    int m_id;
    string m_name;
    float m_dist;
    vector<TreeNode*> m_children;
    TreeNode* m_parent;
};


class Tree
{
public:
    Tree() :
        m_root(NULL)
    {}
    
    ~Tree()
    {
        for (unsigned int i=0; i<m_nodes.size(); i++)
            delete m_nodes[i];
    }
    
    // input/output
    bool Read(FILE *infile);
    TreeNode *ReadNode(FILE *infile, TreeNode *parent, int &depth);
    void Write(FILE *outfile);
    void WriteNode(FILE *outfile, TreeNode *node, int indent);
    void Display(FILE *outfile=stdout, int width=80, float scale=-1);
    
    int *GetSizes();
    
    // accessors
    TreeNode *AddNode(TreeNode* node, TreeNode *parent);    
    inline TreeNode *GetRoot() { return m_root; }
    inline TreeNode *GetNode(int id) { return m_nodes[id]; }
    inline int Size() { return m_nodes.size(); }

private:

    TreeNode *m_root;
    vector<TreeNode*> m_nodes;
};



};
