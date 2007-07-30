/***************************************************************************
* Summon
* Matt Rasmussen
* HashTable.h
*
***************************************************************************/

#ifndef SPIDIR_HASH_TABLE_H
#define SPIDIR_HASH_TABLE_H

#include <list>
#include <string>
#include <vector>


namespace spidir
{

using namespace std;


struct HashInt {
    static unsigned int hash(const int &n) { return (unsigned) n; }
};

struct HashCharStar {
    static unsigned int hash(const char *s)
    {
        unsigned int h = 0, g;
        
        for (; *s; s++) {
            h = (h << 4) + *s;
            if (g = h & 0xF0000000)
                h ^= g >> 24;
            h &= ~g;
        }
        
        return h;
    }    
};

struct HashString {
    static unsigned int hash(string s) 
    {
        return HashCharStar::hash(s.c_str());
    }
};


struct HashPointer {
    static unsigned int hash(const void *p) { return (unsigned int) p; }
};



template <class KeyType, class ValueType, class HashFunc>
class HashTable
{
public:
    HashTable(int size = 100, ValueType null = ValueType()) :
        m_table(NULL),
        m_null(null)
    {
        if (size > 0)
            create(size);
    }
    
    virtual ~HashTable()
    {
        if (m_table)
            delete [] m_table;
    }
    
    void create(int size)
    {
        m_size = size;
        m_table = new ListType [size];
    }
    
    inline ValueType &insert(const KeyType &key, const ValueType &object)
    {
        unsigned int hash = HashFunc::hash(key);
        NodeType node(key, object)
        m_table[hash % m_size].push_back(node);
        return node.second;
    }
    
    void remove(const KeyType &key)
    {
        unsigned int hash = HashFunc::hash(key);
        ListType *lst  = &m_table[hash % m_size];
        
        for (typename ListType::iterator i = lst->begin(); 
             i != lst->end(); i++)
        {
            if ((*i).first == key) {
                lst->erase(i);
                return;
            }
        }
    }
    
    ValueType get(const KeyType &key) {
        unsigned int hash = HashFunc::hash(key);
        ListType *lst = &m_table[hash % m_size];
        
        for (typename ListType::iterator i = lst->begin(); i != lst->end(); i++) {
            if ((*i).first == key) {
                return (*i).second;
            }
        }
        
        return m_null;
    }
    
    bool hasKey(const KeyType &key)
    {
        unsigned int hash = HashFunc::hash(key);
        ListType *lst = &m_table[hash % m_size];
        
        for (typename ListType::iterator i = lst->begin(); i != lst->end(); i++) {
            if ((*i).first == key) {
                return true;
            }
        }
        
        return false;
    }
    
    void keys(vector<KeyType> *keys)
    {
        keys->clear();
        for (int l=0; l<m_size; l++) {
            ListType *lst = &m_table[l];
        
            for (typename ListType::iterator i = lst->begin(); 
                 i != lst->end(); i++) 
            {
                keys->push_back((*i).first);
            }
        }
    }
    
    
    ValueType &operator[](const KeyType &key) {
        unsigned int hash = HashFunc::hash(key);
        ListType *lst = &m_table[hash % m_size];
        
        for (typename ListType::iterator i = lst->begin(); i != lst->end(); i++) {
            if ((*i).first == key) {
                return (*i).second;
            }
        }
        
        // existing key not found, insert new value
        return insert(key, m_null);
    }
    
protected:
    typedef pair<KeyType, ValueType> NodeType;
    typedef list< NodeType > ListType;
    ListType *m_table;
    int m_size;
    ValueType m_null;
};



} // namespace spidir

#endif // SPIDIR_HASH_TABLE_H
