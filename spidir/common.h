#ifndef SPIDIR_COMMON_H
#define SPIDIR_COMMON_H

/*=============================================================================

    SPIDIR - SPecies Informed DIstanced-based Reconstruction
    
    Matt Rasmussen
    Wed Jun 13 22:09:24 EDT 2007


=============================================================================*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string>
#include <vector>

#include "ExtendArray.h"

using namespace std;


// constants
#define PI 3.1415926
#ifndef INFINITY
#   define INFINITY 1e1000
#endif

// indexing a matrix stored as a single array
#define matind(m, i, j) ((m)*(i) + (j))


//=============================================================================
// Math

// computes the log(normalPdf(x | u, s^2))
float normallog(float x, float u, float s);
double gammln(double xx);
float gammalog(float x, float a, float b);
void invertPerm(int *perm, int *inv, int size);

template <class T>
void permute(T* array, int *perm, int size)
{
    ExtendArray<T> tmp(0, size);
    tmp.extend(array, size);
    
    // transfer permutation to temp array
    for (int i=0; i<size; i++)
        tmp[i] = array[perm[i]];
    
    // copy permutation back to original array
    for (int i=0; i<size; i++)
        array[i] = tmp[i];
}


inline float frand()
{
    return rand() / float(RAND_MAX);
}


template <class T>
int findval(T *array, int size, const T &val)
{
    for (int i=0; i<size; i++)
        if (array[i] == val)
            return i;
    return -1;
}



//=============================================================================
// input/output


class BufferedReader
{
public:
    BufferedReader(FILE *stream=NULL, bool autoclose=true) :
        m_stream(stream),
        m_line(0, 10000),
        m_autoclose(autoclose)
    {}
    
    virtual ~BufferedReader()
    {
        if (m_autoclose && m_stream)
            fclose(m_stream);
    }
    
    
    bool open(const char *filename, char *mode, 
              const char *errmsg="cannot read file '%s'\n")
    {
        m_stream = fopen(filename, mode);
        
        if (!m_stream) {
            fprintf(stderr, errmsg, filename);        
            return false;
        }
        return true;
    }
    
    
    char *readLine()
    {
        while (!feof(m_stream)) {
            int pos = m_line.size();
            char *ret = fgets(&(m_line.get()[pos]), 
                              m_line.capacity()-m_line.size(), m_stream);
            int readsize = strlen(&(m_line.get()[pos]));
            
            if (ret == NULL)
                return NULL;
            
            if (m_line.size() + readsize < m_line.capacity() - 1)
                return m_line.get();

            assert(m_line.increaseCapacity());
            m_line.setSize(m_line.size() + readsize);
        }
        
        return NULL;
    }
    
    
    void close()
    {
        if (m_stream)
            fclose(m_stream);
        m_stream = NULL;
    }
    
    
protected:
    FILE *m_stream;
    ExtendArray<char> m_line;
    bool m_autoclose;
};




void printIntArray(int *array, int size);
void printFloatArray(float *array, int size);
bool inChars(char c, const char *chars);
bool chomp(char *str);
vector<string> split(const char *str, const char *delim, bool multiDelim = true);
void printError(const char *fmt, ...);
char readChar(FILE *stream, int &depth);
char readUntil(FILE *stream, string &token, char *stops, int &depth);
string trim(const char *word);



#endif // SPIDIR_COMMON_H
