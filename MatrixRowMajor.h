/* 
 * File:   MatrixRowMajor.h
 * Author: mdesana
 *
 * Created on November 18, 2016, 10:44 AM
 */

#pragma once

#include <assert.h>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <vector>
#include <cstring>

/*a matrix that can be efficiently processed in raster order*/
template <typename T>
class Matrix_rowMajor
{
    //    std::vector< std::vector<T> > m_elements;
    std::vector< std::vector<T> > m_elements;
    int m_rows, m_cols, m_size;

public:

    Matrix_rowMajor()
    {
        m_rows = 0;
        m_cols = 0;
        m_size = 0;
        m_elements = std::vector< std::vector<T> >();
    }

    inline int Size()
    {
        return m_size;
    }

    Matrix_rowMajor(int rows, int cols)
    {
        //        todo redo more efficient
        assert(rows > 0 && cols > 0);
        m_elements = std::vector < std::vector<T> >(rows);
        for (int r = 0; r < rows; r++)
            m_elements[r] = std::vector < T >(cols);
        m_cols = cols;
        m_rows = rows;
        m_size = cols*rows;
    }

    //    Matrix_rowMajor<T>& operator=(const Matrix_rowMajor<T>& rhs)
    //    {
    //        m_size = rhs.m_size;
    //        m_cols = rhs.m_cols;
    //        m_rows = rhs.m_rows;
    //        m_elements = rhs.m_elements;
    //        
    //        return *this;
    //    }

    Matrix_rowMajor GetRowSubset(int startIndex, int stopIndex) const 
    {
        if(stopIndex>nRows())
            stopIndex = nRows();
        int N = stopIndex - startIndex;
        assert(N > 0);
        assert(startIndex >= 0);

        Matrix_rowMajor m(N,nCols());
        int k=0;
        for (int r = startIndex; r < stopIndex; r++)
        {
            const std::vector< T >& row_r = m_elements[r];
            std::vector< T >& other_row_r = m.m_elements[k];
            for (int c = 0; c < m_cols; c++)
                other_row_r[c] = row_r[c];
            k++;
        }
        
        return m;
    }

    const std::vector< T >& GetRow(int r) const
    {
        return m_elements[r];
    }

    std::vector< T >& EditRow(int r)
    {
        return m_elements[r];
    }

    void SetVal(int r, int c, T val)
    {
        m_elements[r][c] = val;
    }

    void Increment(int r, int c, T val)
    {
        m_elements[r][c] += val;
    }

    bool empty() const
    {
        return m_elements.empty();
    }

    T GetVal(int r, int c) const
    {
        //whenever possible, do not use this but GetRow and iterate over elements (it's faster)
        return m_elements[r][c];
    }

    T& EditVal(int r, int c)
    {
        //whenever possible, do not use this but GetRow and iterate over elements (it's faster)
        return m_elements[r][c];
    }

    void ApplyLog()
    {
        for (int r = 0; r < m_rows; r++)
        {
            std::vector< T >& row_r = m_elements[r];
            for (int c = 0; c < m_cols; c++)
                row_r[c] = log(row_r[c]);
        }
    }

    void Sum(const Matrix_rowMajor<T>& other)
    {
        assert(other.nCols() == this->nCols() && other.nRows() == this->nRows());
        for (int r = 0; r < m_rows; r++)
        {
            std::vector< T >& row_r = m_elements[r];
            const std::vector< T >& other_row_r = other.m_elements[r];
            for (int c = 0; c < m_cols; c++)
                row_r[c] += other_row_r[c];
        }
    }

    std::string ToString() const
    {
        std::stringstream s;
        for (int i = 0; i < m_rows; i++)
        {
            s << "  [";
            for (int j = 0; j < m_cols; j++)
            {
                s << GetVal(i, j) << ", ";
            }
            s << "]\n";
        }
        return s.str();
    }

    std::string ToStringExponentiated() const
    {
        std::stringstream s;
        for (int i = 0; i < m_rows; i++)
        {
            s << "  [";
            for (int j = 0; j < m_cols; j++)
            {
                s << exp(GetVal(i, j)) << ", ";
            }
            s << "]\n";
        }
        return s.str();
    }

    Matrix_rowMajor<T> Transposed() const
    {
        Matrix_rowMajor<T> tr(this->nCols(), this->nRows());
        for (int r = 0; r < nRows(); r++)
        {
            for (int c = 0; c < nCols(); c++)
                tr.SetVal(c, r, GetVal(r, c));
        }
        return tr;
    }

    int nRows() const
    {
        return m_rows;
    }

    int nCols() const
    {
        return m_cols;
    }

    void SetVal(T val)
    {
        //whenever possible, do not use this but EditRow and iterate over row elements (it's faster)
        for (int r = 0; r < m_rows; r++)
        {
            std::vector< T >& row_r = m_elements[r];
            for (int c = 0; c < m_cols; c++)
                row_r[c] = val;
        }
    }

    static void Test()
    {
        Matrix_rowMajor<double> mat(2, 3);
        mat.SetVal(0.4);
        std::cout << mat.ToString() << std::endl;
        for (int i = 0; i < mat.nRows(); i++)
            for (int j = 0; j < mat.nCols(); j++)
            {
                //                std::cout << "elem for index " << i << " " << j << " is " << mat.GetIndex(i, j) << std::endl;
                mat.SetVal(i, j, i + j);
            }
        std::cout << mat.ToString() << std::endl;

        Matrix_rowMajor<double> m2;
        m2 = mat;
        m2.SetVal(0, 0, 10);
        std::cout << " values after duplication\n";
        std::cout << mat.ToString() << std::endl;
        std::cout << m2.ToString() << std::endl;
    }
};
//
//
///*a matrix that can be efficiently processed in raster order*/
//template <typename T>
//class Matrix_rowMajor
//{
//    //    std::vector< std::vector<T> > m_elements;
//    std::vector<T> m_elements;
//    int m_rows, m_cols, m_size;
//
//    inline int GetIndex(int r, int c) const
//    {
//        assert(c >= 0 && c < m_cols && r >= 0 && r < m_rows);
//        int ind = r * m_cols + c;
//        assert(ind<m_size);
//        return ind;
//    }
//public:
//
//    Matrix_rowMajor()
//    {
//        m_rows = 0;
//        m_cols = 0;
//        m_size = 0;
//        m_elements = std::vector<T>();
//    }
//
//    inline int Size()
//    {
//        return m_size;
//    }
//
//    Matrix_rowMajor(int rows, int cols)
//    {
//        //        todo redo more efficient
//        assert(rows > 0 && cols > 0);
//        m_elements = std::vector<T>(rows * cols);
//        m_cols = cols;
//        m_rows = rows;
//        m_size = cols*rows;
//    }
//
////    Matrix_rowMajor<T>& operator=(const Matrix_rowMajor<T>& rhs)
////    {
////        m_size = rhs.m_size;
////        m_cols = rhs.m_cols;
////        m_rows = rhs.m_rows;
////        m_elements = rhs.m_elements;
////        
////        return *this;
////    }
//
//    void SetVal(int r, int c, T val)
//    {
//        m_elements[GetIndex(r, c)] = val;
//    }
//
//    void Increment(int r, int c, T val)
//    {
//        m_elements[GetIndex(r, c)] += val;
//    }
//
//    T GetVal(int r, int c) const
//    {
//        //TODO
//        return m_elements[GetIndex(r, c)];
//    }
//
//    void ApplyLog()
//    {
//        for (int i = 0; i < m_size; i++)
//            m_elements[i] = log(m_elements[i]);
//    }
//
//    void Sum(const Matrix_rowMajor<T>& other)
//    {
//        assert(other.nCols() == this->nCols() && other.nRows() == this->nRows());
//        const std::vector<T>& otherElems = other.m_elements;
//        for (int i = 0; i < m_size; i++)
//            m_elements[i] += otherElems[i];
//    }
//
//    std::string ToString() const
//    {
//        std::stringstream s;
//        for (int i = 0; i < m_rows; i++)
//        {
//            s << "  [";
//            for (int j = 0; j < m_cols; j++)
//            {
//                s << GetVal(i, j) << ", ";
//            }
//            s << "]\n";
//        }
//        return s.str();
//    }
//
//    std::vector<T> RowToVector(int r) const
//    {
//        int rowSt = GetIndex(r, 0);
//        std::vector<T> out(m_cols);
//        for (int i = 0; i < m_cols; i++)
//            out[i] = m_elements[i+rowSt];
//        return out;
//    }
//
//    std::string ToStringExponentiated() const
//    {
//        std::stringstream s;
//        for (int i = 0; i < m_rows; i++)
//        {
//            s << "  [";
//            for (int j = 0; j < m_cols; j++)
//            {
//                s << exp(GetVal(i, j)) << ", ";
//            }
//            s << "]\n";
//        }
//        return s.str();
//    }
//
//    Matrix_rowMajor<T> Transposed()
//    {
//        Matrix_rowMajor<T> tr(this->nCols(), this->nRows());
//        for (int r = 0; r < nRows(); r++)
//        {
//            for (int c = 0; c < nCols(); c++)
//                tr.SetVal(c, r, GetVal(r, c));
//        }
//        return tr;
//    }
//
//    int nRows() const
//    {
//        return m_rows;
//    }
//
//    int nCols() const
//    {
//        return m_cols;
//    }
//
//    void SetVal(T val)
//    {
//        //        std::fill(m_elements, m_elements + m_size , val);
//        for (int i = 0; i < m_size; i++)
//            m_elements[i] = val;
//    }
//
//    static void Test()
//    {
//        Matrix_rowMajor<double> mat(2, 3);
//        mat.SetVal(0.4);
//        std::cout << mat.ToString() << std::endl;
//        for (int i = 0; i < mat.nRows(); i++)
//            for (int j = 0; j < mat.nCols(); j++)
//            {
//                std::cout << "elem for index " << i << " " << j << " is " << mat.GetIndex(i, j) << std::endl;
//                mat.SetVal(i, j, i + j);
//            }
//        std::cout << mat.ToString() << std::endl;
//        
//        Matrix_rowMajor<double> m2;
//        m2 = mat;
//        m2.SetVal(0, 0, 10);
//        std::cout<<" values after duplication\n";
//        std::cout << mat.ToString() << std::endl;
//        std::cout << m2.ToString() << std::endl;
//    }
//};