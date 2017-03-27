/* 
 * File:   globals.h
 * Author: mdesana
 *
 * Created on May 8, 2014, 2:33 PM
 */

#pragma once

#include "Params.h" //must be included as first!
#include <assert.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <math.h>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>
#include "Log.h" 
#include <fstream>
#include <list>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <set>
#include "MatrixRowMajor.h"


typedef double Real;

#define RAND_INIT_WEIGHT_MIN 0.25
#define RAND_INIT_WEIGHT_MAX 0.75
#define ROOT_VERTEX 999999999
#define EMPTY_VERTEX 999999998
#define ROOT_VAR -1
#define VAR_NON_INITIALIZED -2
#define NAN_DOUBLE std::numeric_limits<Real>::quiet_NaN()
#define ZERO_LOGVAL -std::numeric_limits<Real>::infinity()

#include "MatrixContiguous.h"
//#define NDEBUG  //if you define NDEBUG, asserts are not evaluated

class SPGMnode;
//We use the property vertex_rank to encode the pointer to the JunctionGraphNode data. Not so nice but works. TODO define a new property with a proper name.


typedef boost::property<boost::vertex_rank_t, SPGMnode*> VertexProperty;
typedef boost::property < boost::edge_weight_t, Real > EdgeWeightProperty;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, VertexProperty, EdgeWeightProperty> Digraph;
typedef boost::graph_traits<Digraph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Digraph>::edge_descriptor Edge;
typedef Digraph::out_edge_iterator out_edge_iterator;
typedef Digraph::in_edge_iterator in_edge_iterator;




//typedef boost::numeric::ublas::matrix<Real> RealMatrix;
//typedef boost::numeric::ublas::matrix<int> IntMatrix;

Real SPN_Rand(Real fMin = 0.f, Real fMax = 1.f);

Real SPN_RandSeed();

Real Sigmoid(Real x);

void WriteToFile(const std::string& str, const std::string& filename);
// 01 Knapsack with int costs. Returns the maximum value that can be put in a knapsack of capacity maxCost and the configuration for the maximum.

Real KnapSack(int maxCost, std::vector<int>& costs, std::vector<Real>& values, std::vector<bool>& chosenItems);

Real AddLog(Real l1, Real l2);

void String2File(const std::string& s,const std::string filename);
std::string GetDateTimeString();

template<typename T>
inline std::string VecToString(const std::vector<T>& v)
{
    std::stringstream s;
    s << "  [";
    for (int i = 0; i < v.size(); i++)
    {
        s << v[i] << ", ";
    }
    s << "]";
    return s.str();
}

template<typename T>
inline std::string VecToStringExponentiated(const std::vector<T>& v)
{
    std::stringstream s;
    s << "  [";
    for (int i = 0; i < v.size(); i++)
    {
        s << exp(v[i]) << ", ";
    }
    s << "]";
    return s.str();
}


inline Real SafeLog(Real d)
{
    if (d == 0)
    {
        return std::numeric_limits<Real>::min();
    }
    else
    {
        return log(d);
    }
}

std::string SetToString(const std::set<int>& mySet);

std::string SetToString(const std::set<Vertex>& mySet);

typedef MatrixContiguous<Real> RealMatrix;

typedef unsigned char byte;

typedef Matrix_rowMajor<byte> ByteMatrix;

class DataMatrix
{
    ByteMatrix m_byRow;
    ByteMatrix m_byCol;
public:

    void InitFromMatbyRow(const ByteMatrix &mat)
    {
        m_byRow = mat;
        m_byCol = mat.Transposed();
    }

    const ByteMatrix& ByCol() const
    {
        return m_byCol;
    }

    ByteMatrix& EditByRow()
    {
        return m_byRow;
    }

    const ByteMatrix& ByRow() const
    {
        return m_byRow;
    }

    int nVars() const
    {
        return m_byRow.nCols();
    }

    int nData() const
    {
        return m_byRow.nRows();
    }

    void ReadFromCSV(std::string filename);

    DataMatrix DataSubset(int startIndex, int stopIndex) const;
};

struct Params
{
    int K;
    int nInsertions; //number of extra edges added to spgm in learning
    Real lambda;
    int nEmExternalIters;
    int nEmInternalIters;
    int EM_iters1;
    int EM_iters2;
    bool verbose;
    bool reproducePaperResults;
    int nIt_EM;
    int nIt_params;
    int nIt_struct;
    bool runAll;
    int NtriesBeforeEarlyStopping;
    std::string filename;
    std::string dataFolder;
    bool em_W;
    bool em_theta;
    Real diricheletPrior;
    Real uniformWeightPrior;
    long int maxBatchMemSize;
    bool fixedSizePerSPGM;
    bool doParamEmInStructLearn;
    int maxNAddedPaths; //for performances. Do not even consider adding more than these paths.  
    int takeNbestMi;
    bool plotTestLL;

    enum InitParams
    {
        init_rand, init_kmeans
    };

    InitParams init;

    Params();
    
    std::string GetUniqueString() const;//gets a short string that should be very likely unique for each setting of parameters 
    
    void PrintHelp();

    void LoadFromArgs(int argc, char** argv);
    std::string ToString() const;
};
