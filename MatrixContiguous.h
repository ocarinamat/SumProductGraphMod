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

#include "Common.h"

/*a matrix that can be efficiently processed in raster order and holds things in a contiguous memory spape (a vector)*/
template <typename T>
class MatrixContiguous
{
    //    std::vector< std::vector<T> > m_elements;
    std::vector<T> m_elements;
    int m_rows, m_cols, m_size;

    inline int ComputeIndex(int r, int c) const
    {
        assert(c >= 0 && c < m_cols && r >= 0 && r < m_rows);
        int ind = r * m_cols + c;
        assert(ind < m_size);
        return ind;
    }
public:

    MatrixContiguous()
    {
        m_rows = 0;
        m_cols = 0;
        m_size = 0;
        m_elements = std::vector<T>();
        
//        std::cout<<"creating empty MatrixContiguous "<<this<<std::endl<<ToString()<<std::endl;
    }
    
//    ~MatrixContiguous()
//    {
//        std::cout<<"deleting MatrixContiguous "<<this<<std::endl<<ToString()<<std::endl;
//    }
    
    inline int Size()
    {
        return m_size;
    }

    MatrixContiguous(int rows, int cols)
    {
        //        todo redo more efficient
        assert(rows > 0 && cols > 0);
        m_elements = std::vector<T>(rows * cols);
        m_cols = cols;
        m_rows = rows;
        m_size = cols*rows;
        
//        std::cout<<"creating  MatrixContiguous"<<m_rows<<" x "<< m_cols<<std::endl;
    }

    void Increment(int r, int c, T val)
    {
        m_elements[ComputeIndex(r, c)] += val;
    }
    
    //    MatrixContiguous<T>& operator=(const MatrixContiguous<T>& rhs)
    //    {
    //        m_size = rhs.m_size;
    //        m_cols = rhs.m_cols;
    //        m_rows = rhs.m_rows;
    //        m_elements = rhs.m_elements;
    //        
    //        return *this;
    //    }

    std::vector<T> ComputeRow(int r) const
    {
        std::vector<T> out(m_cols);
        for (int c = 0; c < m_cols; c++)
            out[c] = m_elements[ComputeIndex(r, c)];
        return out;
    }

    std::vector<T> ComputeColumn(int c) const
    {
        std::vector<T> out(m_rows);
        for (int r = 0; r < m_rows; r++)
            out[r] = m_elements[ComputeIndex(r, c)];
        return out;
    }

    void Sum(const MatrixContiguous<T>& other)
    {
        assert(other.nCols() == this->nCols() && other.nRows() == this->nRows());
        for (int i = 0; i < m_size; i++)
            m_elements[i] += other.m_elements[i];
    }

    void SetVal(int r, int c, T val)
    {
        m_elements[ComputeIndex(r, c)] = val;
    }

    //    void Increment(int r, int c, T val)
    //    {
    //        m_elements[ComputeIndex(r, c)] += val;
    //    }

    T GetVal(int r, int c) const
    {
        //TODO
        return m_elements[ComputeIndex(r, c)];
    }

    T& EditVal(int r, int c)
    {
        //TODO
        return m_elements[ComputeIndex(r, c)];
    }

    bool empty() const
    {
        return m_elements.empty();
    }

    void AddLogVal(int r, int c, const T& logVal)
    {
        int ind = ComputeIndex(r, c);
        T el = m_elements[ind];
        if (el == ZERO_LOGVAL)
            m_elements[ind] = logVal;
        else
            m_elements[ind] = AddLog(el, logVal);
    }

    void ApplyLog()
    {
        for (int i = 0; i < m_size; i++)
            m_elements[i] = log(m_elements[i]);
    }

    bool HasNan() const 
    {
        for (int i = 0; i < m_size; i++)
            if (isnan(m_elements[i]))
                return true;
        return false;
    }

    //    void Sum(const MatrixContiguous<T>& other)
    //    {
    //        assert(other.nCols() == this->nCols() && other.nRows() == this->nRows());
    //        const std::vector<T>& otherElems = other.m_elements;
    //        for (int i = 0; i < m_size; i++)
    //            m_elements[i] += otherElems[i];
    //    }

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

    MatrixContiguous<T> Transposed()
    {
        MatrixContiguous<T> tr(this->nCols(), this->nRows());
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
        //        std::fill(m_elements, m_elements + m_size , val);

        for (int i = 0; i < m_size; i++)
            m_elements[i] = val;
    }

    static void Test()
    {
        MatrixContiguous<double> mat(2, 3);
        mat.SetVal(0.4);
        std::cout << mat.ToString() << std::endl;
        for (int i = 0; i < mat.nRows(); i++)
            for (int j = 0; j < mat.nCols(); j++)
            {
                std::cout << "elem for index " << i << " " << j << " is " << mat.ComputeIndex(i, j) << std::endl;
                mat.SetVal(i, j, i + j);
            }
        std::cout << mat.ToString() << std::endl;

        MatrixContiguous<double> m2;
        m2 = mat;
        m2.SetVal(0, 0, 10);
        std::cout << " values after duplication\n";
        std::cout << mat.ToString() << std::endl;
        std::cout << m2.ToString() << std::endl;
    }
};