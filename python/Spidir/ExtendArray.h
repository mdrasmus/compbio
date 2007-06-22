#ifndef SPIDIR_EXTEND_ARRAY_H
#define SPIDIR_EXTEND_ARRAY_H

#include <assert.h>
#include <algorithm>



/*=============================================================================
    equivalent of realloc for C++
*/
template <class T>
T *resize(T *array, size_t oldsize, size_t newsize)
{
    T *tmp = new T[newsize];
    
    // failed to get memory
    if (!tmp)
        return NULL;
    
    // equivalent to new
    if (oldsize == 0)
        return tmp;
    
    copy(array, array + oldsize, tmp);
   
    delete [] array;
    return tmp;
}


/*=============================================================================
    easy, detachable wrapper for arrays

    len      -- the length of the populated data
    datasize -- the size of the allocated array
*/
template <class ValueType>
class ExtendArray
{
public:
    typedef ValueType* ValuePtrType;

    ExtendArray(int _len=0, int size=0, ValueType *_data=NULL, 
                int _minsize=40) :
        data(_data),
        len(_len),
        datasize(size),
        minsize(_minsize)
    {
        // makesure capacity is atleast length of data
        if (datasize < len)
            datasize = len;
        
        // makesure min allocate size is atleast capacity
        // useful for multiple detachment
        if (minsize < datasize)
            minsize = datasize;
        
        // if no pointer is given to manage, allocate our own
        if (data == NULL && datasize != 0) {
            data = new ValueType [datasize];
        }
    }
    
    ~ExtendArray()
    {
        if (data)
            delete [] data;
    }
    
    ValueType *detach()
    {
        ValueType *ret = data;
        data = NULL;
        len = 0;
        datasize = 0;
        
        return ret;
    }
    
    //=========================================================================
    // capacity management
    
    bool setCapacity(int newsize)
    {
        int oldsize = datasize;
        ValueType *ret = resize(data, oldsize, newsize);
        
        // failed to alloc memory
        if (!ret)
            return false;
        
        data = ret;
        datasize = newsize;
        return true;
    }
    
    bool increaseCapacity()
    {
        int newsize = datasize;
        if (newsize < minsize)
            newsize = minsize;
        newsize *= 2;
        return setCapacity(newsize);
    }
    
    bool ensureSize(int needed)
    {
        if (needed <= datasize)
            return true;
        
        int newsize = needed;
        if (newsize < minsize)
            newsize = minsize;
        while (newsize < needed)
            newsize *= 2;
        
        return setCapacity(newsize);
    }
    
    inline int capacity()
    {
        return datasize;
    }
    
    
    //=========================================================================
    // data access
    void append(const ValueType &val)
    {
        assert(ensureSize(len + 1));
        data[len++] = val;
    }
    
    void extend(ValueType *vals, int nvals)
    {
        assert(ensureSize(len + nvals));
        
        for (int i=0; i<nvals; i++)
            data[len++] = vals[i];
    }
    
    ValueType pop()
    {
        assert(len > 0);
        return data[--len];
    }
    
    inline ValueType &operator[](const int i)
    {
        return data[i];
    }
    
    inline ValueType *get()
    {
        return data;
    }
    
    // easy access to underling data
    operator ValuePtrType()
    {
        return data;
    }
    
    
    inline int size()
    {
        return len;
    }
    
    inline void setSize(int _size)
    {
        len = _size;
    }
    
protected:    
    ValueType *data;
    int len;
    int datasize;
    int minsize;
};



#endif // SPIDIR_EXTEND_ARRAY_H
