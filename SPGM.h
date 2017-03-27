/* 
 * File:   SPGM.h
 * Author: mdesana
 *
 * Created on October 17, 2016, 3:17 PM
 */

#pragma once

#include "Common.h"
#include "Factor.h"
#include "SPGMnodes.h"

class SPGM
{
public:
    bool m_verbose;
    bool m_debug;
    int m_maxMemForBatch;

    //Real LearnChowLiuSPGM(Dataset data);
    //Real LearnMixMaxSpanningTreesSPGM(Dataset data);

    SPGM() : m_verbose(false), m_debug(false), m_dataDim(2), m_mem_PerSample_perEval(-1)
    {
        m_maxMemForBatch = std::numeric_limits<int>::infinity();
    }
    ~SPGM(); //TODO
    SPGM(const SPGM& other);
    SPGM& operator= ( const SPGM & other);  
    void ComputeScope();
    void ComputeVParents();

    const std::set<int>& GetVParents(Vertex v) const
    {
        std::map<Vertex, std::set<int> >::const_iterator it = m_vParentsMap.find(v);
        assert(it != m_vParentsMap.end() && "item not found");
        return it->second;
    }

    const std::set<Vertex>& GetVParentVertices(Vertex v) const
    {
        std::map<Vertex, std::set<Vertex> >::const_iterator it = m_vParentsVerticesMap.find(v);
        assert(it != m_vParentsVerticesMap.end() && "item not found");
        return it->second;
    }

    Real MeanLL(const DataMatrix& data)
    {
        const std::vector<Real>& ll = Eval(data);
        Real sum = 0.0;
        for (int i = 0; i < ll.size(); i++)
            sum += ll[i];
        return sum / ll.size();
    }

    std::vector<Real> Eval(const DataMatrix& data);
    
    Real EvalAndDerivate_batch();

    SPGM EM(const DataMatrix& training, const DataMatrix& validation, int maxIters = 200, int NtriesBeforeEarlyStopping = 3, Real uniformWeightPrior = 0, bool weightUpdate = true, bool factorsUpdate = true);

    Vertex AddVNode(int var)
    {
        Vertex newV = NewVertex();
        Vnode* vn = new Vnode(var, newV);
        boost::property_map<Digraph, boost::vertex_rank_t>::type nodeData = get(boost::vertex_rank, m_graph);
        nodeData[newV] = vn;
        return newV;
    }

    Vertex AddProdNode()
    {
        Vertex newV = NewVertex();
        ProdNode* pn = new ProdNode(newV);
        boost::property_map<Digraph, boost::vertex_rank_t>::type nodeData = get(boost::vertex_rank, m_graph);
        nodeData[newV] = pn;
        return newV;
    }

    Vertex AddSumNode()
    {
        Vertex newV = NewVertex();
        SumNode* sn = new SumNode(newV);
        boost::property_map<Digraph, boost::vertex_rank_t>::type nodeData = get(boost::vertex_rank, m_graph);
        nodeData[newV] = sn;
        return newV;
    }

    void AddEdge(Vertex parentId, Vertex childId)
    {
        //todo use weighted edges
        SPGMnode* n = EditNode(parentId);
        if ((n->GetType() == SPGMnode::vnode) && (GetChildrenVertices(parentId).size() > 0))
        {
            std::cerr << "Vnodes must have at most 1 children";
        }

        boost::add_edge(parentId, childId, m_graph);
        m_finalized = false;
    }

    //    void AddSumEdge(Vertex parentId, Vertex childId, Real weight)
    //    {
    //        //todo use weighted edges
    //        SPGMnode* n = EditNode(parentId);
    //        assert(n->GetType() == SPGMnode::sum); //only sum nodes have weighted edges
    //        assert(weight>=0);
    //        boost::add_edge(parentId, childId, weight, m_graph);
    //        m_finalized = false;
    //    }

    void RemoveEdge(Vertex parentId, Vertex childId)
    {
        boost::remove_edge(parentId, childId, m_graph);
        m_finalized = false;
    }

    std::list<Vertex> GetChildrenVertices(Vertex v) const;
    std::list<Vertex> GetParentVertices(Vertex v) const;

    void SetFactor(Vertex child, Vertex parent, const CPTfactor& f);
    void SetFactor2Root(Vertex child, const CPTfactor& f);

    SPGMnode* EditNode(Vertex v);

    const SPGMnode* GetNodeInOrder(int i) const
    {
        return m_invTopOrder[i];
    }

    void WriteGraphViz(std::string filename);

    void ClearMessages();
    void Finalize();

    static void Test();

    std::string ToString(bool printParams);

    std::vector<const SPGMnode*> GetInvTopOrder() const
    {
        assert(m_finalized);
        std::vector<const SPGMnode*> out(m_invTopOrder.size());
        for (int i = 0; i < m_invTopOrder.size(); i++)
            out[i] = m_invTopOrder[i];
        return out;
    };

    void ComputeTopOrder();

    int GetNMessages() const
    {
        return m_nMsgs;
    }

    void EM_resetBetas();

    int GetNNodes() const
    {
        return m_graph.m_vertices.size();
    }

    void ComputeDerivatives(const DataMatrix& data, Real rootLogDerivative, const std::vector<Real>& externalSampleLogP);

    void ComputeDerivatives(const DataMatrix& data)
    {
        const RealMatrix& rootMat = m_invTopOrder[m_invTopOrder.size() - 1]->EditMsgLogVals(ROOT_VERTEX);
        const std::vector<Real>& sampleLogP = rootMat.ComputeColumn(0);
        ComputeDerivatives(data, 0.0, sampleLogP);
    }

    void EM_updateStep(Real uniformWeightPrior, bool weightUpdate, bool factorsUpdate);

    void StructLearn_MixMaxSpanTrees();

    const std::list<CPTfactor>& GetFactors()
    {
        return m_factors;
    }
    
    int GetMemPerSample() const
    {
        assert(m_finalized);
        return m_mem_PerSample_perEval;
    }

private:
    //    std::map<Vertex, std::list<int> > m_vertex2VparentsMap; //vars in the vParents for each node
    std::list<CPTfactor> m_factors;
    std::vector<SPGMnode*> m_invTopOrder;
    std::map<Vertex, std::set<int> > m_scopeMap;
    std::map<Vertex, std::set<int> > m_vParentsMap; //set of vParent variables for each vertex
    std::map<Vertex, std::set<Vertex> > m_vParentsVerticesMap; //vertices corresponding to each vParent
    std::set<int> m_allVars;
    int m_nMsgs; 
    long m_mem_PerSample_perEval;

    //todo make const
    int m_dataDim; //note, for now we only accept variables of the same size. (and of size 2...))

    int m_oldBatchSize;
    Digraph m_graph;

    void DeepCopyAux(const SPGM& other);
    
    void InstantiateMessages(int N); //set the number of samples in the data to be processed
    
    void ComputeNMessages();

    Vertex NewVertex()
    {
        m_finalized = false;
        //add a vertex to the graph
        Vertex u;
        u = boost::add_vertex(m_graph);
        return u;
    }

    bool m_finalized;

};

