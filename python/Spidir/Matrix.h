/*******************************************************************************
 SINDER - Matt Rasmussen
 Matrix.h
 11/19/05

*******************************************************************************/

#include <stdio.h>




template <class ValueType>
class Matrix
{
public:
    Matrix(int nrows, int ncols) :
        m_nrows(nrows),
        m_ncols(ncols)
    {
        m_data = new (ValueType*)[nrows];
        for (int i=0; i<nrows; i++)
            m_data[i] = new ValueType[ncols];
    }
    
    ~Matrix()
    {
        for (int i=0; i<m_nrows; i++)
            delete [] m_data[i];
        delete [] m_data;
    }
    
    void SetAll(ValueType val)
    {
        for (int i=0; i<m_nrows; i++)
            for (int j=0; j<m_ncols; j++)
                m_data[i][j] = val;
    }
    
    inline int NumRows() { return m_nrows; }
    inline int NumCols() { return m_ncols; }
    
    inline ValueType* operator[](const int i)
    {
        return m_data[i];
    }
    
    inline ValueType **getMatrix()
    {
        return m_data;
    }
    
private:
    int m_nrows;
    int m_ncols;
    ValueType **m_data;
};


