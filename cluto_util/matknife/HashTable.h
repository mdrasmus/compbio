/***************************************************************************
* Vistools
* Matt Rasmussen
* HashTable.h
*
***************************************************************************/

#ifndef HASH_TABLE_H
#define HASH_TABLE_H

#include <list>
#include <string>
#include <vector>


using namespace std;


struct HashInt {
    static unsigned int Hash(const int &n) { return (unsigned) n; }
};

struct HashCharStar {
    static unsigned int Hash(unsigned char *str)
    {
        unsigned int hash = 5381;
        int c;

        while ((c = *str++))
            hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

        return hash;
    }
};

struct HashString {
    static unsigned int Hash(string s) 
    {
        return HashCharStar::Hash((unsigned char*) s.c_str());
    }
};



template <class KeyType, class ValueType, class HashFunc>
class HashTable
{
public:
    HashTable(int size = 20, ValueType null = ValueType()) :
        m_table(NULL),
        m_null(null)
    {
        if (size > 0)
            Create(size);
    }
    
    virtual ~HashTable()
    {
        if (m_table)
            delete [] m_table;
    }
    
    void Create(int size)
    {
        m_size = size;
        m_table = new ListType [size];
    }
    
    int Size()
    {
        int size = 0;
        for (int i=0; i<m_size; i++)
            size += m_table[i].size();
        return size;
    }
    
    inline ValueType &Insert(const KeyType &key, const ValueType &object)
    {
        unsigned int hash = HashFunc::Hash(key) % m_size;
        m_table[hash].push_back(Node(key, object));
        return m_table[hash].back().value;
    }
    
    void Remove(const KeyType &key)
    {
        unsigned int hash = HashFunc::Hash(key);
        ListType *lst  = &m_table[hash % m_size];
        
        for (typename ListType::iterator i = lst->begin(); 
             i != lst->end(); i++)
        {
            if ((*i).key == key) {
                lst->erase(i);
                return;
            }
        }
    }
    
    ValueType &Get(const KeyType &key) {
        unsigned int hash = HashFunc::Hash(key);
        ListType *lst = &m_table[hash % m_size];
        
        for (typename ListType::iterator i = lst->begin(); i != lst->end(); i++) {
            if ((*i).key == key) {
                return (*i).value;
            }
        }
        
        return Insert(key, m_null);
    }
    
    bool HasKey(const KeyType &key)
    {
        unsigned int hash = HashFunc::Hash(key);
        ListType *lst = &m_table[hash % m_size];
        
        for (typename ListType::iterator i = lst->begin(); i != lst->end(); i++) {
            if ((*i).key == key) {
                return true;
            }
        }
        
        return false;
    }
    
    void Keys(vector<KeyType> *keys)
    {
        keys->clear();
        for (int l=0; l<m_size; l++) {
            ListType *lst = &m_table[l];
        
            for (typename ListType::iterator i = lst->begin(); 
                 i != lst->end(); i++) 
            {
                keys->push_back((*i).key);
            }
        }
    }
    
    inline ValueType &operator[](const KeyType &key)
    {
        return Get(key);
    }
    
    struct Node {
        Node(KeyType _key, ValueType _val) :
            key(_key), value(_val)
        {}
        KeyType key;
        ValueType value;
    };
    
    typedef list< Node > ListType;
    
    class Iterator {
    public:
        Iterator(HashTable *_hash) :
            hash(_hash),
            arrayPos(0),
            listPos(_hash->m_table[0].begin())
        {
            SetNext();
        }
        
        Node Next()
        {
            if (HasMore()) {
                Node *node = &(*listPos);
                SetNext();
                return *node;
            }
            return Node(KeyType(), hash->m_null);
        }
        
        inline bool HasMore()
        {
            return (arrayPos < hash->m_size);
        }
        
        void Reset()
        {
            arrayPos = 0;
            listPos = hash->m_table[0]->begin();
            SetNext();
        }       
        
    protected:
        inline void NextNonEmptyList()
        {
            // find next non-empty list
            do {
                arrayPos++;
            } while(arrayPos < hash->m_size && 
                    hash->m_table[arrayPos].size() == 0);
            if (arrayPos < hash->m_size)
                listPos = hash->m_table[arrayPos].begin();
        }
        
        void SetNext()
        {
            if (listPos == hash->m_table[arrayPos].end()) {
                NextNonEmptyList();
            } else {
                listPos++;
                if (listPos == hash->m_table[arrayPos].end())
                    NextNonEmptyList();
            }
        }
        
        friend class HashTable;
        HashTable *hash;
        int arrayPos;
        typename ListType::iterator listPos;
    };
    
    Iterator Begin()
    {
        return Iterator(this);
    }
    
    
protected:
    
    ListType *m_table;
    int m_size;
    ValueType m_null;
};


#endif
