
#include "SPGM.h"
#include "opencv2/opencv.hpp"
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <time.h>

using namespace boost;
using namespace std;

void SPGM::SetFactor(Vertex child, Vertex parent, const CPTfactor& f)
{
    m_finalized = false;
    SPGMnode* n = EditNode(child);
    assert(n->GetType() == SPGMnode::vnode);
    SPGMnode* p = EditNode(parent);
    assert(p->GetType() == SPGMnode::vnode);
    Vnode* vChild = (Vnode*) n;
    Vnode* vParent = (Vnode*) p;
    int paVar = vParent->GetVar();

    //if there is not a message to parent variable, create it
    Message* m = vChild->GetMessage(parent);
    if (m == NULL)
    {
        vChild->AddMessage(parent, paVar);
        if (m_verbose)
            std::cout << "  Adding a message from vertex " << child << " to vertex " << parent << std::endl;
    }
    else
    {
        if (m_verbose)
            std::cerr << "message from vertex " << child << " to vertex " << parent << " already there\n";
    }

    //set the factor
    m_factors.push_back(f);
    CPTfactor* fPointer = &m_factors.back();
    vChild->SetFactor(fPointer, parent);
    //    bool found = false;
    //    for (list<CPTfactor*>::iterator it = m_factors.begin(); it != m_factors.end(); it++)
    //    {
    //        if ((*it) == f)
    //        {
    //            found = true;
    //            break;
    //        }
    //    }
    //    if (!found)
    //        m_factors.push_back(f);
}

void SPGM::DeepCopyAux(const SPGM& other)
{
    //copy the non dynamic part
    m_allVars = other.m_allVars;
    m_dataDim = other.m_dataDim;
    m_debug = other.m_debug;
    m_factors = other.m_factors;
    m_finalized = other.m_finalized;
    m_graph = other.m_graph;
    m_oldBatchSize = -1;
    m_scopeMap = other.m_scopeMap;
    m_vParentsMap = other.m_vParentsMap;
    m_vParentsVerticesMap = other.m_vParentsVerticesMap;
    m_verbose = other.m_verbose;
    m_nMsgs = other.m_nMsgs;
    m_mem_PerSample_perEval = other.m_mem_PerSample_perEval;

    //now do the dynamic part
    map<const CPTfactor*, CPTfactor*> fMap;
    list<CPTfactor>::iterator itFthis = m_factors.begin();
    list<CPTfactor>::const_iterator itFother = other.m_factors.begin();
    for (; itFthis != m_factors.end(); itFthis++, itFother++)
    {
        //record the address of old factor corresponding to address of new factor. 
        CPTfactor const* otherFacPointer = &(*itFother);
        fMap[otherFacPointer] = &(*itFthis);
    }
    //    map<const SPGMnode*, SPGMnode*> nMap;
    m_invTopOrder = std::vector<SPGMnode*>(other.GetNNodes());
    for (int i = 0; i < other.GetNNodes(); i++)
    {
        const SPGMnode* nOther = other.GetNodeInOrder(i);
        if (nOther->GetType() == SPGMnode::vnode)
        {
            Vnode* vn = (Vnode*) nOther;
            m_invTopOrder[i] = vn->CreateCopy(fMap);
            //            nMap[nOther] = vn;
        }
        else if (nOther->GetType() == SPGMnode::sum)
        {
            SumNode* vn = (SumNode*) nOther;
            m_invTopOrder[i] = vn->CreateCopy();
            //            nMap[nOther] = vn;
        }
        else
        {
            ProdNode* vn = (ProdNode*) nOther;
            m_invTopOrder[i] = vn->CreateCopy();
            //            nMap[nOther] = vn;
        }
    }

    //now change the SPGMnode pointer information in the graph data
    boost::property_map<Digraph, boost::vertex_rank_t>::type nodeData = get(boost::vertex_rank, m_graph);
    for (int i = 0; i < m_invTopOrder.size(); i++)
    {
        SPGMnode* currNode = m_invTopOrder[i];
        //        Vertex p = m_invTopOrder[i]->GetVertex();
        //        SPGMnode* otherNode = nodeData[p];
        nodeData[currNode->GetVertex()] = currNode;
    }

    //now add the children messages

    //todo this part should only be called once we actually evaluate the SPGM 
    for (int i = 0; i < m_invTopOrder.size(); i++)
    {
        SPGMnode* currNode = m_invTopOrder[i];
        Vertex p = m_invTopOrder[i]->GetVertex();
        //        SPGMnode* currNode = nodeData[p];
        currNode->ClearChildren();
        const set<Vertex>& vPaVertices = m_vParentsVerticesMap[currNode->GetVertex()];

        //#######################    now add children messages 
        out_edge_iterator begin, end;
        boost::tie(begin, end) = boost::out_edges(p, m_graph);
        for (out_edge_iterator currEdge = begin; currEdge != end; ++currEdge)
        {
            SPGMnode* child = nodeData[boost::target(*currEdge, m_graph)];
            if (currNode->GetType() == SPGMnode::vnode)
            {
                //find message from child to currNode
                Message* m = child->GetMessage(currNode->GetVertex());
                if (m == NULL)
                    printf("ERROR: CHILD %s has no message to %s\n", child->ToString().c_str(), currNode->ToString().c_str());

                //pass it to EVERY output message (including root)
                for (set<Vertex>::const_iterator it = vPaVertices.begin(); it != vPaVertices.end(); it++)
                {
                    Vertex currPa = *it;
                    currNode->AddChildMsg(m, currPa);
                }

                if (m_verbose)
                    printf("  ADD MSG FROM CHILD %s to %s\n", child->ToString().c_str(), currNode->ToString().c_str());
            }
            else
            {
                set<Vertex>& vPaVertices = m_vParentsVerticesMap[currNode->GetVertex()];
//                for (set<Vertex>::const_iterator it = vPaVertices.begin(); it != vPaVertices.end(); it++)
//                {
//                    //find all messages from child to the vParent vertices of currNode
//                    Vertex vPa = *it;
//                    SPGMnode* vpaNode = EditNode(vPa);
//                    Message* m = child->GetMessage(vPa);
//                    if (m == NULL)
//                        printf("ERROR: CHILD %s has no message to %s\n", child->ToString().c_str(), vpaNode->ToString().c_str());
//
//                    currNode->AddChildMsg(m, vPa);
//                    if (m_verbose)
//                        printf("  ADD MSG FROM CHILD %s to %s passing through node %s. \n", child->ToString().c_str(), vpaNode->ToString().c_str(), currNode->ToString().c_str());
//                }
                
                for (set<Vertex>::const_iterator it = vPaVertices.begin(); it != vPaVertices.end(); it++)
                {
                    //find all messages from child to the vParent vertices of currNode
                    Vertex vPa = *it;
                    Message* m = child->GetMessage(vPa);
                    if (vPa != ROOT_VERTEX)
                    {
                        SPGMnode* vpaNode = EditNode(vPa);
                        if (m == NULL)
                            printf("ERROR: CHILD %s has no message to %s\n", child->ToString().c_str(), vpaNode->ToString().c_str());
                        if (m_verbose)
                            printf("  ADD MSG FROM CHILD %s to %s passing through node %s. \n", child->ToString().c_str(), vpaNode->ToString().c_str(), currNode->ToString().c_str());
                    }
                    else
                    {
                        if (m == NULL)
                            printf("ERROR: CHILD %s has no message to ROOT\n", child->ToString().c_str());
                        if (m_verbose)
                            printf("  ADD MSG FROM CHILD %s to ROOT passing through node %s. \n", child->ToString().c_str(), currNode->ToString().c_str());
                    }
                    currNode->AddChildMsg(m, vPa);
                }
                
            }
        }

        if (begin != end)
            assert(currNode->NChildren() > 0);
        if (currNode->GetType() != SPGMnode::vnode)
            assert(currNode->NChildren() > 0);
    }
}

SPGM& SPGM::operator=(const SPGM & other)
{
    //first destroy dynamic existing things
    for (int i = 0; i < m_invTopOrder.size(); i++)
        delete m_invTopOrder[i];

    DeepCopyAux(other);
    return *this;
}

SPGM::SPGM(const SPGM& other) : m_verbose(false), m_debug(false), m_dataDim(2)
{
    //initializes the SPGM performing a deep copy of the other.
    //evaluated messages are not included in the copy
    DeepCopyAux(other);
}

std::list<Vertex> SPGM::GetChildrenVertices(Vertex v) const
{
    list<Vertex> children;
    out_edge_iterator begin, end;
    boost::tie(begin, end) = boost::out_edges(v, m_graph);
    for (out_edge_iterator currEdge = begin; currEdge != end; ++currEdge)
    {
        Vertex c = boost::target(*currEdge, m_graph);
        children.push_back(c);
    }
    return children;
}

std::list<Vertex> SPGM::GetParentVertices(Vertex v) const
{
    list<Vertex> parents;
    in_edge_iterator begin, end;
    boost::tie(begin, end) = boost::in_edges(v, m_graph);
    for (in_edge_iterator currEdge = begin; currEdge != end; ++currEdge)
    {
        Vertex c = boost::source(*currEdge, m_graph);
        parents.push_back(c);
    }
    return parents;
}

void SPGM::SetFactor2Root(Vertex child, const CPTfactor& f)
{
    m_finalized = false;
    SPGMnode* n = EditNode(child);
    assert(n->GetType() == SPGMnode::vnode);

    Vnode* vChild = (Vnode*) n;
    Vertex paVertex = ROOT_VERTEX;

    //if there is not a message to parent variable, create it
    Message* m = vChild->GetMessage(paVertex);
    if (m == NULL)
    {
        vChild->AddMessage(paVertex, ROOT_VAR);
        if (m_verbose)
            std::cout << "  Adding a message from vertex " << child << " to root " << std::endl;
    }
    else
    {
        if (m_verbose)
            std::cerr << "message from vertex " << child << " to root already there\n";
    }

    //set the factor
    m_factors.push_back(f);
    CPTfactor* fPointer = &m_factors.back();
    vChild->SetFactor(fPointer, paVertex);

    //    vChild->SetFactor(f, paVertex);
    //    bool found = false;
    //    for (list<CPTfactor*>::iterator it = m_factors.begin(); it != m_factors.end(); it++)
    //    {
    //        if ((*it) == f)
    //        {
    //            found = true;
    //            break;
    //        }
    //    }
    //    if (!found)
    //        m_factors.push_back(f);
}

void SPGM::ComputeTopOrder()
{
    //compute the topological order
    if (m_verbose)
        cout << "Computing inverse topological order" << endl;
    boost::property_map<Digraph, boost::vertex_rank_t>::type nodeData = get(boost::vertex_rank, m_graph);
    std::vector< Vertex > c;
    boost::topological_sort(m_graph, std::back_inserter(c));
    m_invTopOrder = std::vector<SPGMnode*>(c.size());
    for (int i = 0; i < c.size(); i++)
    {
        Vertex& p = c[i];
        SPGMnode* currNode = nodeData[p];
        if (m_verbose)
            printf("  INV TOP ORDER %d: %s\n", i, currNode->ToString().c_str());
        m_invTopOrder[i] = currNode;
    }
}

void SPGM::ClearMessages()
{
    //clear the children nodes (will be reassigned later based on the graph)) and the sum and prod messages
    boost::property_map<Digraph, boost::vertex_rank_t>::type nodeData = get(boost::vertex_rank, m_graph);
    for (int i = 0; i < m_invTopOrder.size(); i++)
    {
        Vertex p = m_invTopOrder[i]->GetVertex();
        SPGMnode* currNode = nodeData[p];
        currNode->ClearChildren();

        if ((currNode->GetType() == SPGMnode::sum) || (currNode->GetType() == SPGMnode::prod))
            currNode->ClearMessages();
    }
}

void SPGM::Finalize()
{
    if (m_finalized)
        return; //do nothing

    ComputeTopOrder();
    //compute vparents 
    ComputeVParents();
    boost::property_map<Digraph, boost::vertex_rank_t>::type nodeData = get(boost::vertex_rank, m_graph);

    //clear the children nodes (will be reassigned later based on the graph)) and the sum and prod messages
    ClearMessages();

    // collect all vars in the SPGM.
    m_allVars = set<int>();
    for (int i = 0; i < m_invTopOrder.size(); i++)
    {
        Vertex p = m_invTopOrder[i]->GetVertex();
        SPGMnode* currNode = nodeData[p];
        if (currNode->GetType() == SPGMnode::vnode)
        {
            Vnode* vn = (Vnode*) currNode;
            m_allVars.insert(vn->GetVar());
        }
    }

    //add the children to each node, based on the graph edges. Also add the children messages. 
    if (m_verbose)
        cout << "Adding children messages to all SPGMnodes, and output messages to prod and sum nodes." << endl;
    for (int i = 0; i < m_invTopOrder.size(); i++)
    {
        Vertex p = m_invTopOrder[i]->GetVertex();
        SPGMnode* currNode = nodeData[p];
        currNode->ClearChildren();
        const set<Vertex>& vPaVertices = m_vParentsVerticesMap[currNode->GetVertex()];

        //##################   create messages for prod and sum nodes to every vParent vertex
        int currVertex = currNode->GetVertex();
        if (currNode->GetType() == SPGMnode::sum || currNode->GetType() == SPGMnode::prod)
        {
            assert(currNode->Nmessages() == 0); //should not have been added previously

            for (set<Vertex>::const_iterator it = vPaVertices.begin(); it != vPaVertices.end(); it++)
            {
                Vertex currPa = *it;
                if (currPa == ROOT_VERTEX)
                {
                    currNode->AddMessage(ROOT_VERTEX, ROOT_VAR);
                    if (m_verbose)
                        cout << "  Adding message from vertex " << currVertex << " to root " << endl;
                }
                else
                {
                    assert(EditNode(currPa)->GetType() == SPGMnode::vnode);
                    Vnode* vn = (Vnode*) EditNode(currPa);
                    int nStatesTo = vn->GetMessage(vn->GetOutputMessageVertices().front())->m_factor->NIn(); //TODO better system to know var size
                    currNode->AddMessage(currPa, vn->GetVar());
                    if (m_verbose)
                        cout << "  Adding message from vertex " << currVertex << " to vertex " << currPa << endl;
                }
            }
        }

        //#######################    now add children messages 
        out_edge_iterator begin, end;
        boost::tie(begin, end) = boost::out_edges(p, m_graph);
        for (out_edge_iterator currEdge = begin; currEdge != end; ++currEdge)
        {
            SPGMnode* child = nodeData[boost::target(*currEdge, m_graph)];
            if (currNode->GetType() == SPGMnode::vnode)
            {
                //find message from child to currNode
                Message* m = child->GetMessage(currNode->GetVertex());
                if (m == NULL)
                    printf("ERROR: CHILD %s has no message to %s\n", child->ToString().c_str(), currNode->ToString().c_str());

                //pass it to EVERY output message (including root)
                for (set<Vertex>::const_iterator it = vPaVertices.begin(); it != vPaVertices.end(); it++)
                {
                    Vertex currPa = *it;
                    currNode->AddChildMsg(m, currPa);
                }

                if (m_verbose)
                    printf("  ADD MSG FROM CHILD %s to %s\n", child->ToString().c_str(), currNode->ToString().c_str());
            }
            else
            {
                set<Vertex>& vPaVertices = m_vParentsVerticesMap[currNode->GetVertex()];
                for (set<Vertex>::const_iterator it = vPaVertices.begin(); it != vPaVertices.end(); it++)
                {
                    //find all messages from child to the vParent vertices of currNode
                    Vertex vPa = *it;
                    Message* m = child->GetMessage(vPa);
                    if (vPa != ROOT_VERTEX)
                    {
                        SPGMnode* vpaNode = EditNode(vPa);
                        if (m == NULL)
                            printf("ERROR: CHILD %s has no message to %s\n", child->ToString().c_str(), vpaNode->ToString().c_str());
                        if (m_verbose)
                            printf("  ADD MSG FROM CHILD %s to %s passing through node %s. \n", child->ToString().c_str(), vpaNode->ToString().c_str(), currNode->ToString().c_str());
                    }
                    else
                    {
                        if (m == NULL)
                            printf("ERROR: CHILD %s has no message to ROOT\n", child->ToString().c_str());
                        if (m_verbose)
                            printf("  ADD MSG FROM CHILD %s to ROOT passing through node %s. \n", child->ToString().c_str(), currNode->ToString().c_str());
                    }
                    currNode->AddChildMsg(m, vPa);
                }
            }
        }
    }
    m_finalized = true;
    ComputeNMessages();
    //    ComputeScope();
}

void SPGM::ComputeVParents()
{
    m_vParentsMap = map<Vertex, std::set<int> > ();
    m_vParentsVerticesMap = map<Vertex, std::set<Vertex> >();

    if (m_verbose)
        cout << "ComputeVParents:\n";

    for (int i = m_invTopOrder.size() - 1; i >= 0; i--)
    {
        SPGMnode* node = m_invTopOrder[i];
        m_vParentsMap[node->GetVertex()] = set<int>();
        m_vParentsVerticesMap[node->GetVertex()] = set<Vertex>();
    }

    //set the root 
    Vertex topVertex = m_invTopOrder[m_invTopOrder.size() - 1]->GetVertex();
    m_vParentsMap[topVertex].insert(ROOT_VAR);
    m_vParentsVerticesMap[topVertex].insert(ROOT_VERTEX);

    for (int i = m_invTopOrder.size() - 1; i >= 0; i--)
    {
        SPGMnode* node = m_invTopOrder[i];
        Vertex currVertex = node->GetVertex();

        const list<Vertex>& children = GetChildrenVertices(currVertex);
        const set<int>& vPas = m_vParentsMap[currVertex];
        const set<Vertex>& vPaVertices = m_vParentsVerticesMap[currVertex];

        if (node->GetType() == SPGMnode::vnode)
        {
            //add to the child the currVertex to vPas
            Vnode* vn = (Vnode*) node;
            assert(children.size() <= 1);
            if (!children.empty())
            {
                set<int>& chVpas = m_vParentsMap[children.front()];
                chVpas.insert(vn->GetVar());
                m_vParentsVerticesMap[children.front()].insert(node->GetVertex());
                if (chVpas.size() != 1)
                    cerr << "ERROR: child vParents should have 1 element only: " << SetToString(chVpas) << endl;
            }
        }
        else
        {
            //pass the vParents to the children 
            for (list<Vertex>::const_iterator chIt = children.begin(); chIt != children.end(); chIt++)
            {
                set<int>& chVpas = m_vParentsMap[*chIt];
                for (set<int>::const_iterator it = vPas.begin(); it != vPas.end(); it++)
                    chVpas.insert(*it);
                set<Vertex>& chvPaVertices = m_vParentsVerticesMap[*chIt];
                for (set<Vertex>::const_iterator it = vPaVertices.begin(); it != vPaVertices.end(); it++)
                    chvPaVertices.insert(*it);
            }
        }

        if (m_verbose)
        {
            cout << "  vParents of node " << node->ToString() << ": " << SetToString(vPas) << "\n     vPa vertices: " << SetToString(vPaVertices) << endl;
        }

    }
}

void SPGM::ComputeScope()
{
    assert(m_finalized);

    //compute and check the scope of each node, in inverse topological order
    m_scopeMap = map<Vertex, set<int> > ();

    for (int i = m_invTopOrder.size() - 1; i >= 0; i--)
    {
        SPGMnode* node = m_invTopOrder[i];
        m_scopeMap[node->GetVertex()] = set<int>();
    }


    if (m_verbose)
        cout << "CheckScope:\n";
    for (int i = 0; i < m_invTopOrder.size(); i++)
    {
        SPGMnode* node = m_invTopOrder[i];
        m_scopeMap[node->GetVertex()] = set<int>();
        set<int>& scope = m_scopeMap[node->GetVertex()];
        const list<Vertex>& children = GetChildrenVertices(node->GetVertex());

        if (m_verbose)
        {
            cout << "  scope of node " << node->ToString() << ": ";
        }

        if (node->GetType() == SPGMnode::vnode)
        {
            //scope is scope of child plus internal var
            //internal var must not be in child
            Vnode* vn = (Vnode*) m_invTopOrder[i];
            assert(children.size() <= 1);
            if (!children.empty())
                scope = m_scopeMap[children.front()];
            int var = vn->GetVar();
            if (scope.find(var) != scope.end())
                cerr << "\nERROR: scope " << SetToString(scope) << " contains " << var << endl;
            scope.insert(var);
        }
        else if (node->GetType() == SPGMnode::sum)
        {
            //scope is the scope of the first child
            //scopes must be the same
            scope = m_scopeMap[children.front()];
            list<Vertex>::const_iterator chIt = children.begin();
            chIt++; //start from second child
            for (; chIt != children.end(); chIt++)
            {
                const set<int>& childScope = m_scopeMap[*chIt];
                if (childScope.size() != scope.size())
                {
                    cerr << "\nERROR: child scope " << SetToString(childScope) << " different from scope " << SetToString(scope) << endl;
                    exit(1);
                }
                for (set<int>::iterator it = childScope.begin(); it != childScope.end(); it++)
                {
                    if (scope.find(*it) == scope.end())
                        cerr << "\nERROR: child scope " << SetToString(scope) << " does not contain " << *it << endl;
                }
            }
        }
        else if (node->GetType() == SPGMnode::prod)
        {
            //scope is the union of child scopes
            //child scopes must be disjoint
            scope = m_scopeMap[children.front()];
            list<Vertex>::const_iterator chIt = children.begin();
            chIt++; //start from second child
            for (; chIt != children.end(); chIt++)
            {
                const set<int>& childScope = m_scopeMap[*chIt];
                for (set<int>::iterator it = childScope.begin(); it != childScope.end(); it++)
                {
                    if (scope.find(*it) != scope.end())
                        cerr << "\nERROR: child scope " << SetToString(scope) << " contains " << *it << endl;
                    scope.insert(*it);
                }
            }
        }

        if (m_verbose)
        {
            cout << SetToString(scope) << endl;
        }

    }
}

SPGMnode* SPGM::EditNode(Vertex v)
{
    boost::property_map<Digraph, boost::vertex_rank_t>::type nodeData = get(boost::vertex_rank, m_graph);
    return nodeData[v];
}

SPGM::~SPGM()
{
    //    for (std::list<CPTfactor*>::iterator it = m_factors.begin(); it != m_factors.end(); it++)
    //        delete *it;
    //    if (m_verbose)
    //        cout<<"deleting SPGM "<<this<<endl;
    for (int i = 0; i < m_invTopOrder.size(); i++)
        delete m_invTopOrder[i];
}

void SPGM::Test()
{
    //create a spgm A + (  B   *   (C1 -> D) + (C2 -> D) )
    SPGM* spgmP = new SPGM();
    SPGM& spgm = *spgmP;
    spgm.m_verbose = false;

    int dim = 2;
    CPTfactor f_a(1, dim);
    CPTfactor f_ba(dim, dim);
    CPTfactor f_ca1(dim, dim);
    CPTfactor f_ca2(dim, dim);
    CPTfactor f_dc1(dim, dim);
    CPTfactor f_dc2(dim, dim);
    //    
    Vertex a = spgm.AddVNode(0);
    //    Vertex prod = spgm.AddProdNode();
    Vertex b = spgm.AddVNode(1);
    //        Vertex sum = spgm.AddSumNode();
    //        Vertex c1 = spgm.AddVNode(2);
    //        Vertex c2 = spgm.AddVNode(2);
    //        Vertex d = spgm.AddVNode(3);

    spgm.SetFactor2Root(a, f_a);
    spgm.SetFactor(b, a, f_ba);
    //        spgm.SetFactor(c1, a, f_ca1);
    //        spgm.SetFactor(c2, a, f_ca2);
    //        spgm.SetFactor(d, c1, f_dc1);
    //        spgm.SetFactor(d, c2, f_dc2);
    //

    spgm.AddEdge(a, b);

    //    spgm.AddEdge(a, prod);
    //    spgm.AddEdge(prod, b);
    //    spgm.AddEdge(prod, sum);
    //    spgm.AddEdge(sum, c1);
    //    spgm.AddEdge(sum, c2);
    //    spgm.AddEdge(c1, d);
    //    spgm.AddEdge(c2, d);
    //
    spgm.WriteGraphViz("SPGMtest");
    //
    spgm.Finalize();
    //
    //    SumNode* sn = (SumNode*) spgm.EditNode(sum);
    //    sn->SetWeightsRandom();
    //
    ByteMatrix data(1, 2);
    data.SetVal(0, 0, 0);
    data.SetVal(0, 1, 1);
    //    data.SetVal(0, 2, 1);
    //    data.SetVal(0, 3, 0);
    DataMatrix dataMat;
    dataMat.InitFromMatbyRow(data);
    const vector<Real>& llVec = spgm.Eval(dataMat);
    cout << "Probs " << VecToStringExponentiated(llVec) << endl;

    //    spgm.EM(data, data, 200, 3, 0, true, true);

    SPGM spgmCopy(spgm);
    double ll = spgm.MeanLL(dataMat);
    double llCopy = spgmCopy.MeanLL(dataMat);

    cout << "ORIGINAL:\n" << spgm.ToString(true);
    cout << "COPY:\n" << spgmCopy.ToString(true);
    printf("\n#### RESULTS: SPGM LL %.12f SPGMcopy LL %.12f difference %.12f\n", ll, llCopy, ll - llCopy);

    spgm.m_factors.front().SetRandomWeights();
    ll = spgm.MeanLL(dataMat);

    delete spgmP;

    llCopy = spgmCopy.MeanLL(dataMat);
    printf("\n#### AFTER CHANGING ORIGINAL FACTORS: SPGM LL %.12f SPGMcopy LL %.12f difference %.12f\n", ll, llCopy, ll - llCopy);

    //    printf("\n#### Checking Memory leaks (memory usage should not increase) \n");
    for (int i = 0; i < 100; i++)
    {
        //        spgmCopy = spgmCopy;
        //        SPGM other(spgmCopy);
        //        cout<<"other LL"<<other.MeanLL(dataMat)<<endl;
        //        spgmCopy = SPGM(spgmCopy);
        //        llCopy = spgmCopy.MeanLL(dataMat);
    }
}

void SPGM::ComputeNMessages()
{
    assert(m_finalized);
    m_nMsgs = 0;
    int theSize = sizeof (Real);
    m_mem_PerSample_perEval = 0;
    for (int i = 0; i < m_invTopOrder.size(); i++)
    {
        SPGMnode* node = m_invTopOrder[i];

        const vector<Message*>& msgs = node->GetMsgs();
        m_nMsgs += msgs.size();
        for (int m = 0; m < msgs.size(); m++)
        {
            if (msgs[m]->m_toVertex == ROOT_VERTEX)
                m_mem_PerSample_perEval += theSize;
            else
                m_mem_PerSample_perEval += theSize*m_dataDim;
        }
    }
    m_mem_PerSample_perEval *= 2;
}

void SPGM::InstantiateMessages(int N)
{
    if (N != m_oldBatchSize)
    {
        for (int i = 0; i < m_invTopOrder.size(); i++)
        {
            SPGMnode* node = m_invTopOrder[i];
            node->CreateMessages(N, m_dataDim);
        }
        m_oldBatchSize = N;
    }
}

std::vector<Real> SPGM::Eval(const DataMatrix & data)
{
    //data must be arranged each sample in column
    assert(data.nVars() == m_allVars.size() && "sample dimension is different to number of variables in the SPGM");
    int N = data.nData();

    //    int nMsgs = this->GetNMessages();
    //    int memRequest = nMsgs * N*m_dataDim; //todo root messages are overestimated, since they are only one row 
    //
    //    if (memRequest > m_maxMemForBatch)
    //        DataMatrix dataBatch;

    InstantiateMessages(N);
    assert(m_finalized);
    for (int i = 0; i < m_invTopOrder.size(); i++)
    {
        SPGMnode* node = m_invTopOrder[i];
        if (node->GetType() == SPGMnode::vnode)
        {
            Vnode* vn = (Vnode*) m_invTopOrder[i];
            vn->Eval(data);
        }
        else
            node->Eval();

    }
    if (m_verbose)
    {
        std::cout << "SPGM_eval " << std::endl;
        for (int i = 0; i < m_invTopOrder.size(); i++)
        {
            SPGMnode* node = m_invTopOrder[i];
            std::cout << node->ToStringMessages(false);
        }
        std::cout << "END SPGM_eval\n\n";
    }
    //return root message
    const RealMatrix& rootMat = m_invTopOrder[m_invTopOrder.size() - 1]->EditMsgLogVals(ROOT_VERTEX);
    return rootMat.ComputeColumn(0);
}

void SPGM::ComputeDerivatives(const DataMatrix& data, Real rootLogDerivative, const std::vector<Real>& sampleLogP)
{
    assert(m_finalized);
    m_invTopOrder[m_invTopOrder.size() - 1]->EditLogDspnDmessage(ROOT_VERTEX).SetVal(rootLogDerivative);

    //pass the derivatives root to leaves
    for (int i = m_invTopOrder.size() - 1; i >= 0; i--)
    {
        SPGMnode* node = m_invTopOrder[i];
        if (node->GetType() == SPGMnode::vnode)
        {
            Vnode* vn = (Vnode*) m_invTopOrder[i];
            vn->PassDerivative_andStackBeta(data, sampleLogP);
        }
        else if (node->GetType() == SPGMnode::sum)
        {
            SumNode* sn = (SumNode*) m_invTopOrder[i];
            sn->PassDerivative_andStackBeta(sampleLogP);
        }
        else if (node->GetType() == SPGMnode::prod)
        {
            node->PassDerivative();
        }

    }
}

SPGM SPGM::EM(const DataMatrix& train, const DataMatrix& valid, int maxIters, int NtriesBeforeEarlyStopping, Real uniformWeightPrior, bool weightUpdate, bool factorsUpdate)
{
    int N = train.nData();
    int worseThanBestValidLLcount = 0;
    clock_t t1, t2;

    if (m_verbose)
        cout << "\n#####learnWithEM\n";

    //compute initial validation LL and save initial params
    if (m_verbose)
        cout << "\n##computing validation LL\n";
    //        SavedParams bestParams = new SavedParams();
    Real bestValidLL = MeanLL(valid);

    if (m_verbose)
        cout << " initial validation LL " << bestValidLL << endl;

    int bestIter = 0;
    SPGM bestSpgm(*this); //save params

    int iter;

    Real superOldTrainLL = 1;

    for (iter = 1; iter < maxIters; iter++)
    {
        if (m_verbose)
            cout << "\n##EM iter " << iter << endl;

        //reset the infos for the E step (beta in sum nodes and factors)
        for (int i = 0; i < m_invTopOrder.size(); i++)
        {
            SPGMnode* node = m_invTopOrder[i];
            if (node->GetType() == SPGMnode::sum)
            {
                SumNode* sn = (SumNode*) node;
                sn->ResetBeta();
            }
        }
        for (list<CPTfactor>::iterator it = m_factors.begin(); it != m_factors.end(); it++)
            (*it).ResetBeta();

        //perform forward and backward pass for all samples in training set, saving the info for the E step
        t1 = clock();
        const vector<Real>& ll = Eval(train);
        t2 = clock();
        float diff1 = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
        t1 = clock();
        ComputeDerivatives(train);
        t2 = clock();
        float diff2 = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;

        Real oldTrainLL = 0;
        for (int n = 0; n < ll.size(); n++)
            oldTrainLL += ll[n];
        oldTrainLL /= ll.size();

        t1 = clock();
        //then perform the M step for sum nodes and factors
        if (weightUpdate)
        {
            for (int i = 0; i < m_invTopOrder.size(); i++)
            {
                SPGMnode* node = m_invTopOrder[i];
                if (node->GetType() == SPGMnode::sum)
                {
                    SumNode* sn = (SumNode*) node;
                    {
                        sn->EMupdate(uniformWeightPrior);
                    }
                }
            }
        }
        if (factorsUpdate)
        {
            for (list<CPTfactor>::iterator it = m_factors.begin(); it != m_factors.end(); it++)
                (*it).EMupdate(uniformWeightPrior);
        }
        t2 = clock();
        float diff3 = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;

        if (m_verbose)
            cout << "\n##computing validation LL\n";

        //evaluate validation LL
        Real validLL = MeanLL(valid);

        if (m_debug && superOldTrainLL != 1 && (superOldTrainLL - oldTrainLL) > 0.00000000000001)
        {
            printf("Warning: training LL decreased from %.12f to %.12f (diff %.13f). This should not happen.\n", superOldTrainLL, oldTrainLL, oldTrainLL - superOldTrainLL);
        }
        superOldTrainLL = oldTrainLL;

        //save best parameters and check convergence
        if (validLL > bestValidLL)
        {
            worseThanBestValidLLcount = 0;
            bestValidLL = validLL;
            bestSpgm = SPGM(*this); //save params
            bestIter = iter;
        }
        else
        {
            worseThanBestValidLLcount++;
            if (worseThanBestValidLLcount > NtriesBeforeEarlyStopping)
            {
                printf("EM terminated at iteration %d by validation LL convergence\n", iter);
                break;
            }
        }

        if (m_verbose)
            cout << "\n######ITER RESULTS\n";
        printf("  EM iter %d: oldTrainLL %f, bestValidLL %f, validLL %f (eval %.1fs, der %.1fs, upd %.1fs)\n", iter, oldTrainLL, bestValidLL, validLL, diff1, diff2, diff3);
        if (m_verbose)
            cout << "\n\n";

    }
    if (iter == maxIters)
    {
        printf("EM terminated at iteration %d by reaching maximum iterations\n", iter - 1);
    }

    //set the best params to the SPN 
    //    (*this) = SPGM(bestSpgm); //save params
    Real currentValidLL = bestSpgm.MeanLL(valid);
    if (bestValidLL != currentValidLL)
    {
        printf("ERROR: mismatch in best parameters LL: bestValidLL %f vs currentValidLL %f\n", bestValidLL, currentValidLL);
    }
    printf("best EM validation LL %f was found at iter %d\n\n", currentValidLL, bestIter);
    return bestSpgm;
}

string SPGM::ToString(bool printParams)
{
    assert(m_finalized);
    stringstream ss;
    ss << "SPGM::ToString(). Nodes in inv topol order: \n";
    for (int i = m_invTopOrder.size() - 1; i >= 0; i--)
    {
        SPGMnode* node = m_invTopOrder[i];
        ss << node->ToString() << endl;
        ss << node->ToStringMessages(printParams);
        if (node->GetType() == SPGMnode::sum)
        {
            SumNode* sn = (SumNode*) node;
            ss << sn->ToStringWeights();
        }
        ss << endl;
    }
    ss << "END OF SPGM::ToString()\n";
    return ss.str();
}

void SPGM::EM_updateStep(Real uniformWeightPrior, bool weightUpdate, bool factorsUpdate)
{
    if (weightUpdate)
    {
        for (int i = 0; i < m_invTopOrder.size(); i++)
        {
            SPGMnode* node = m_invTopOrder[i];
            if (node->GetType() == SPGMnode::sum)
            {
                SumNode* sn = (SumNode*) node;
                {
                    sn->EMupdate(uniformWeightPrior);
                }
            }
        }
    }
    if (factorsUpdate)
    {
        for (list<CPTfactor>::iterator it = m_factors.begin(); it != m_factors.end(); it++)
            (*it).EMupdate(uniformWeightPrior);
    }
}

void SPGM::EM_resetBetas()
{
    for (int i = 0; i < m_invTopOrder.size(); i++)
    {
        SPGMnode* node = m_invTopOrder[i];
        if (node->GetType() == SPGMnode::sum)
        {
            SumNode* sn = (SumNode*) node;
            sn->ResetBeta();
        }
    }
    for (list<CPTfactor>::iterator it = m_factors.begin(); it != m_factors.end(); it++)
        (*it).ResetBeta();
}

/* 
 * writing to a graphViz file
 */

template <class Property>
class NodeWriter
{
public:

    NodeWriter(Property _name) : m_prop(_name)
    {
    }

    template <class VertexOrEdge>
    void operator()(std::ostream& out, const VertexOrEdge& v) const
    {
        SPGMnode* n = get(m_prop, v);
        out << "[fixedsize=true]";

        if (n->GetType() == SPGMnode::prod)
        {
            stringstream ss;
            ss << "x[" << n->GetVertex() << "]";
            out << "[label=" << escape_dot_string(ss.str()) << "]";
            out << "[shape=" << "circle" << "]";
            //            out << "[fillcolor=" << "gray" << "]";
            out << "[width=" << "0.5" << "]";
        }
        else if (n->GetType() == SPGMnode::sum)
        {
            stringstream ss;
            ss << "+[" << n->GetVertex() << "]";
            out << "[label=" << escape_dot_string(ss.str()) << "]";
            out << "[shape=" << "circle" << "]";
            //            out << "[bg=" << "gray" << "]";
            out << "[width=" << "0.5" << "]";
        }
        else if (n->GetType() == SPGMnode::vnode)
        {
            Vnode* vn = (Vnode*) n;
            stringstream ss;
            ss << vn->GetVar() << "[" << vn->GetVertex() << "]";
            out << "[label=" << escape_dot_string(ss.str()) << "]";
            out << "[shape=" << "circle" << "]";
            out << "[width=" << "0.7" << "]";
        }
    }
private:
    Property m_prop;
};

template <class Property>
inline NodeWriter<Property>
MakeNodeWriter(Property n)
{
    return NodeWriter<Property > (n);
}

//the following helper class sets the graph properties in graphviz... 

struct GraphWriter
{

    void operator()(std::ostream&) const
    {
    }

    template <class VorE>
    void operator()(std::ostream& out, const VorE&) const
    {
        out << "splines = line\n";
        out << "ranksep = \"1.1 equally\"\n";
        out << "nodesep=\" 1 \"\n";
    }
};

void SPGM::WriteGraphViz(string name)
{
    property_map<Digraph, vertex_rank_t>::type nodeData = boost::get(vertex_rank, m_graph);
    name += ".dot";
    std::ofstream ofs(name.c_str());
    boost::write_graphviz(ofs, m_graph, MakeNodeWriter(nodeData), GraphWriter());
}
