/* 
 * File:   LearnSPGM.h
 * Author: mdesana
 *
 * Created on November 1, 2016, 8:19 PM
 */

#pragma once

#include "SPGM.h"
#include "SPGM_mixture.h"

struct ProbabilitiesStruct
{
    RealMatrix p00;
    RealMatrix p01;
    RealMatrix p10;
    RealMatrix p11;
    std::vector<Real> p0;
    std::vector<Real> p1;
    std::vector<Real> c0;
    std::vector<Real> c1;
};

struct WeightedEdgeInfo
{
    int a;
    int b;
    Real w;
    bool operator<(const WeightedEdgeInfo& right) const
    {
        return w < right.w;
    }
};

struct NodeInfo
{
    int Vertex;
    int depth;
    const SPGMnode* spgmNode;
    int var;
};

typedef std::map<int, std::list<int> > VarToChildrenVarMap;
typedef std::map<Vertex, std::list<Vertex> > VertexToChildrenVerticesMap;
typedef std::map<int, std::list<Vertex> > VarToChildrenMap;
typedef std::map<int, NodeInfo> varToNodeMap;

struct PathInfo
{
    int rootVar;
    int spgmId;
    std::vector<int> m_vars; //all vars in path
    WeightedEdgeInfo minEdge; //min edge in path
    VarToChildrenVarMap varChdMap;
    VertexToChildrenVerticesMap vertChildMap;
    WeightedEdgeInfo edgeToAdd;
    std::string ToString();
    bool operator<(const PathInfo& right) const;
};

struct MixMax_auxStruct
{
    std::list<PathInfo> pathInfosToAdd;
    std::list<PathInfo> selectedPathInfos;
    Matrix_rowMajor<Real> mi;
    varToNodeMap varToNodeMapInTreeSPGM;
    VertexToChildrenVerticesMap vertexToChdInTreeSPGM;
    ProbabilitiesStruct probs;
};

void ComputeMutualInfo_binary(const DataMatrix& data, std::vector<Real> weights, Matrix_rowMajor<Real>& outputMI, ProbabilitiesStruct& outputProbs, Real diricheletPrior, bool verbose = false);
CPTfactor ComputeCPT(int childVar, int parentVar, const ProbabilitiesStruct& p, bool verbose);
void ChowLiuTree(SPGM& out, const Matrix_rowMajor<Real>& mi, const ProbabilitiesStruct& p, bool verbose = false);
void LearnMixMaxSpanningTreesSPGM(SPGM& outSpgm, const DataMatrix& data, const Params& p);
void LearnMixMaxSpanningTreesSPGM(SPGM& outSpgm, const DataMatrix& data, std::vector<Real> weights, const Params& p);

MixMax_auxStruct FindCandidatePaths(SPGM& outSpgm, const DataMatrix& data, std::vector<Real> weights, const Params& p);
std::list<PathInfo> SelectPaths(std::list<PathInfo>& pathInfosToAdd, const Params& p, bool multiplyPerNmix);
void AddPaths(SPGM& outSpgm, MixMax_auxStruct& auxStruct, const Params& p);

void LearnMixEM(const Params& p, Real* out_trainLL = NULL, Real* out_validLL = NULL, Real* out_testLL = NULL);


void LearnMixMaxSpanningTreesSPGM_factored(SPGM& outSpgm, const DataMatrix& data, std::vector<Real> weights, const Params& p);