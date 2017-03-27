
#include "LearnSPGM.h"
#include <boost/graph/kruskal_min_spanning_tree.hpp>

using namespace std;

using namespace boost;
typedef adjacency_list < vecS, vecS, undirectedS, no_property, property < edge_weight_t, Real > > Ugraph;
typedef graph_traits < Ugraph >::edge_descriptor UEdge;
typedef std::pair<int, int> E;


typedef map<int, list<int> > VarToChildrenVarMap;
typedef map<Vertex, list<Vertex> > VertexToChildrenVerticesMap;
typedef map<int, list<Vertex> > VarToChildrenMap;

string VarToChildrenMapToString(VarToChildrenMap& varChdMap)
{
    stringstream ss;
    for (VarToChildrenMap::iterator it = varChdMap.begin(); it != varChdMap.end(); it++)
    {
        int v = it->first;
        const list<Vertex>& vl = it->second;
        ss << "  {" << v << " : ";
        for (list<Vertex>::const_iterator it2 = vl.begin(); it2 != vl.end(); it2++)
            ss << *it2 << ", ";
        ss << "}\n";
    }
    return ss.str();
}

string VarToChildrenVarMapToString(VarToChildrenVarMap& varChdMap)
{
    stringstream ss;
    for (VarToChildrenVarMap::iterator it = varChdMap.begin(); it != varChdMap.end(); it++)
    {
        int v = it->first;
        const list<int>& vl = it->second;
        ss << "  {" << v << " : ";
        for (list<int>::const_iterator it2 = vl.begin(); it2 != vl.end(); it2++)
            ss << *it2 << ", ";
        ss << "}\n";
    }
    return ss.str();
}

struct NodeInfo
{
    int Vertex;
    int depth;
    const SPGMnode* spgmNode;
    int var;
};

struct WeightedEdgeInfo
{
    int a;
    int b;
    Real w;
};

struct PathInfo
{
    int rootVar;
    vector<int> m_vars; //all vars in path
    WeightedEdgeInfo minEdge; //min edge in path
    VarToChildrenVarMap varChdMap;
    VertexToChildrenVerticesMap vertChildMap;
    WeightedEdgeInfo edgeToAdd;

    string ToString()
    {
        stringstream ss;
        ss << "PATH with root var " << rootVar << " vars: " << VecToString(m_vars) << "\nVertex adj matrix: \n";
        for (VertexToChildrenVerticesMap::iterator it = vertChildMap.begin(); it != vertChildMap.end(); it++)
        {
            Vertex v = it->first;
            const list<Vertex>& vl = it->second;
            ss << "  {" << v << " : ";
            for (list<Vertex>::const_iterator it2 = vl.begin(); it2 != vl.end(); it2++)
                ss << *it2 << ", ";
            ss << "}\n";
        }
        ss << "Var adj matrix: \n";
        for (VarToChildrenVarMap::iterator it = varChdMap.begin(); it != varChdMap.end(); it++)
        {
            int v = it->first;
            const list<int>& vl = it->second;
            ss << "  {" << v << " : ";
            for (list<int>::const_iterator it2 = vl.begin(); it2 != vl.end(); it2++)
                ss << *it2 << ", ";
            ss << "}\n";
        }
        return ss.str();
    }

    bool operator<(const PathInfo& right) const
    {
        return edgeToAdd.w < right.edgeToAdd.w;
    }
};

struct PathCoupleInfo{
    PathInfo original;
    PathInfo new;
};

typedef map<int, NodeInfo> varToNodeMap;

Real SafeLog(Real d)
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

void ComputeMutualInfo_binary(const DataMatrix& data, std::vector<Real> weights, Matrix_rowMajor<Real>& outputMI, ProbabilitiesStruct& outputProbs, Real diricheletPrior, bool verbose)
{

    int Nsamples = data.nData();
    int Nnodes = data.nVars();

    assert(diricheletPrior >= 0);

    if (weights.empty())
    {
        weights = std::vector<Real>(Nsamples);
        for (int i = 0; i < Nsamples; i++)
        {
            weights[i] = 1;
        }
    }
    else
    {
        assert(weights.size() == Nsamples);
        Real sumw = 0;
        for (int i = 0; i < Nsamples; i++)
        {
            sumw += weights[i];
        }
        for (int i = 0; i < Nsamples; i++)
        {
            weights[i] = Nsamples * weights[i] / sumw;
        }
    }

    //    long startTime;
    //    long endTime;
    //    Real duration;

    // COMPUTE COUNTS (WEIGHTED))
    Matrix_rowMajor<Real> c00(Nnodes, Nnodes);
    Matrix_rowMajor<Real> c01(Nnodes, Nnodes);
    Matrix_rowMajor<Real> c10(Nnodes, Nnodes);
    Matrix_rowMajor<Real> c11(Nnodes, Nnodes);
    c00.SetVal(0);
    c01.SetVal(0);
    c10.SetVal(0);
    c11.SetVal(0);
    vector<Real> c0(Nnodes);
    vector<Real> c1(Nnodes);
    for (int v = 0; v < Nnodes; v++)
    {
        c0[v] = 0.0;
        c1[v] = 0.0;
    }

    //    startTime = System.nanoTime();
    for (int i = 0; i < Nsamples; i++)
    {
        Real weight = weights[i];
        //convert to bool so memory representation is compact and more cacheable
        //        vector<bool> sample = vector<bool>(Nnodes);
        //        int dataJ;
        const vector<byte>& dataRow = data.ByRow().GetRow(i);
        //        for (int j = 0; j < Nnodes; j++)
        //        {
        //            dataJ = data.ByRow().GetVal(i, j);
        //
        //            if (dataJ > 0)
        //            {
        //                sample[j] = true;
        //                c1[j] += weight;
        //            }
        //            else
        //            {
        //                sample[j] = false;
        //                c0[j] += weight;
        //            }
        //        }
        //update counts (TODO CHECK BOTTLENECK)
        for (int v = 0; v < Nnodes; v++)
        {
            if (dataRow[v])
            {
                c1[v] += weight;

                for (int w = v + 1; w < Nnodes; w++)
                {
                    if (dataRow[w])
                        c11.Increment(v, w, weight);
                    else
                        c10.Increment(v, w, weight);
                }
            }
            else
            {
                c0[v] += weight;

                for (int w = v + 1; w < Nnodes; w++)
                {
                    if (dataRow[w])
                        c01.Increment(v, w, weight);
                    else
                        c00.Increment(v, w, weight);
                }
            }
        }
    }

    //    endTime = System.nanoTime();
    //    duration = (endTime - startTime) / 1000000000.0;
    //    if (m_printTimes)
    //    {
    //        printf("computing counts time: %f sec\n", duration);
    //    }
    //    startTime = System.nanoTime();


    //COMPUTE EMPIRICAL PROBABILITIES (WEIGHTED) 
    outputProbs.p00 = RealMatrix(Nnodes, Nnodes);
    outputProbs.p01 = RealMatrix(Nnodes, Nnodes);
    outputProbs.p10 = RealMatrix(Nnodes, Nnodes);
    outputProbs.p11 = RealMatrix(Nnodes, Nnodes);
    outputProbs.p0 = vector<Real>(Nnodes);
    outputProbs.p1 = vector<Real>(Nnodes);

    RealMatrix& p00 = outputProbs.p00;
    RealMatrix& p01 = outputProbs.p01;
    RealMatrix& p10 = outputProbs.p10;
    RealMatrix& p11 = outputProbs.p11;
    vector<Real>& p0 = outputProbs.p0;
    vector<Real>& p1 = outputProbs.p1;

    p00.SetVal(0);
    p01.SetVal(0);
    p10.SetVal(0);
    p11.SetVal(0);

    vector<Real> h(Nnodes);
    for (int v = 0; v < Nnodes; v++)
    {
        p0[v] = 0.0;
        p1[v] = 0.0;
        h[v] = NAN;
    }
    for (int v = 0; v < Nnodes; v++)
    {
        //compute probabilities
        p0[v] = (c0[v] + diricheletPrior) / (Nsamples + diricheletPrior * 2);
        p1[v] = (c1[v] + diricheletPrior) / (Nsamples + diricheletPrior * 2);
        h[v] = -p0[v] * SafeLog(p0[v]) - p1[v] * SafeLog(p1[v]);

        for (int w = v + 1; w < Nnodes; w++)
        {
            Real z = c00.GetVal(v, w) + c01.GetVal(v, w) + c10.GetVal(v, w) + c11.GetVal(v, w) + 4 * diricheletPrior;
            p00.SetVal(v, w, (c00.GetVal(v, w) + diricheletPrior) / z);
            p01.SetVal(v, w, (c01.GetVal(v, w) + diricheletPrior) / z);
            p10.SetVal(v, w, (c10.GetVal(v, w) + diricheletPrior) / z);
            p11.SetVal(v, w, (c11.GetVal(v, w) + diricheletPrior) / z);

            //these can be computed simply reversing the items
            p00.SetVal(w, v, p00.GetVal(v, w));
            p01.SetVal(w, v, p10.GetVal(v, w));
            p10.SetVal(w, v, p01.GetVal(v, w));
            p11.SetVal(w, v, p11.GetVal(v, w));

            if (verbose)
            {
                printf("counts for (%d,%d) = [0,0]:%f, [0,1]:%f, [1,0]:%f, [1,1]:%f \n", v, w, c00.GetVal(v, w), c01.GetVal(v, w), c10.GetVal(v, w), c11.GetVal(v, w));
                printf("p for (%d,%d) = p[0,0]:%f, p[0,1]:%f, p[1,0]:%f, p[1,1]:%f \n", v, w, p00.GetVal(v, w), p01.GetVal(v, w), p10.GetVal(v, w), p11.GetVal(v, w));
            }
        }
    }


    //    endTime = System.nanoTime();
    //    duration = (endTime - startTime) / 1000000000.0;
    //    if (m_printTimes)
    //    {
    //        printf("computing probabilities time: %f sec\n", duration);
    //    }
    //    startTime = System.nanoTime();

    //COMPUTE MUTUAL INFO 
    outputMI = Matrix_rowMajor<Real>(Nnodes, Nnodes);
    assert(Nnodes == outputMI.nCols());
    for (int v = 0; v < Nnodes; v++)
    {
        outputMI.SetVal(v, v, h[v]);

        for (int w = v + 1; w < Nnodes; w++)
        {
            Real h_vw = -p00.GetVal(v, w) * SafeLog(p00.GetVal(v, w))
                    - p01.GetVal(v, w) * SafeLog(p01.GetVal(v, w))
                    - p10.GetVal(v, w) * SafeLog(p10.GetVal(v, w))
                    - p11.GetVal(v, w) * SafeLog(p11.GetVal(v, w));
            outputMI.SetVal(v, w, h[v] + h[w] - h_vw);

            if (isnan(outputMI.GetVal(v, w)))
            {
                cerr << "there is a NaN";
            }
            if (verbose)
            {
                printf("Mutual information(%d,%d) = %f\n", v, w, outputMI.GetVal(v, w));
            }
        }
    }
    //complete the symmetric part
    for (int v = 0; v < Nnodes; v++)
    {
        for (int w = 0; w < v; w++)
        {
            outputMI.SetVal(v, w, outputMI.GetVal(w, v));
        }
    }
}

void ChowLiuTree(SPGM& outSpgm, const Matrix_rowMajor<Real>& mi, const ProbabilitiesStruct& p, bool verbose)
{

    //COMPUTE MAXIMUM SPANNING TREE (by computing minimum spanning tree on graph with minus mutual information as edge weight)
    int num_nodes = mi.nCols();
    assert(mi.nCols() == mi.nRows());

    vector<E> edge_array = vector<E>(((num_nodes * num_nodes) - num_nodes) / 2);
    vector<Real> edge_weights = vector<Real>(((num_nodes * num_nodes) - num_nodes) / 2);
    int ind = 0;
    for (int i = 0; i < num_nodes; i++)
    {
        for (int j = i + 1; j < num_nodes; j++)
        {
            edge_array[ind] = E(i, j);
            edge_weights[ind] = -mi.GetVal(i, j);
            //            cout<<"W"<<edge_weights[ind]<<endl;
            ind++;
        }
    }
    assert(ind == edge_array.size() && ind == edge_weights.size());

    Ugraph fullyConnected(edge_array.begin(), edge_array.end(), edge_weights.begin(), num_nodes);

    property_map < Ugraph, edge_weight_t >::type weight = get(edge_weight, fullyConnected);
    std::vector < UEdge > mst;

    kruskal_minimum_spanning_tree(fullyConnected, std::back_inserter(mst));

    if (verbose)
    {
        std::cout << "Print the edges in the MST:" << std::endl;
        for (std::vector < UEdge >::iterator ei = mst.begin();
                ei != mst.end(); ++ei)
        {
            std::cout << source(*ei, fullyConnected) << " <--> " << target(*ei, fullyConnected)
                    << " with weight of " << weight[*ei]
                    << std::endl;
        }

        std::ofstream fout("mst.dot");
        fout << "graph A {\n"
                << " rankdir=LR\n"
                << " size=\"3,3\"\n"
                << " ratio=\"filled\"\n"
                << " edge[style=\"bold\"]\n" << " node[shape=\"circle\"]\n";
        graph_traits<Ugraph>::edge_iterator eiter, eiter_end;
        for (tie(eiter, eiter_end) = edges(fullyConnected); eiter != eiter_end; ++eiter)
        {

            if (std::find(mst.begin(), mst.end(), *eiter) != mst.end())
            {
                fout << source(*eiter, fullyConnected) << " -- " << target(*eiter, fullyConnected);
                fout << "[color=\"black\", label=\"" << get(edge_weight, fullyConnected, *eiter) << "\"];\n";
            }
            //            else{fout << "[color=\"black\", label=\"" << get(edge_weight, mst, *eiter) << "\"];\n";
            //                fout << "[color=\"gray\", label=\"" << get(edge_weight, mst, *eiter)
            //                << "\"];\n";}
        }
        fout << "}\n";
    }


    //create SPGM representing the mst
    outSpgm = SPGM();
    Vertex spgmRootV = outSpgm.AddVNode(0);
    queue<Vertex> toProcess;
    toProcess.push(0);

    //visit nodes in breadth first search and at each node add its children in the SPGM 
    vector<bool> included = vector<bool>(num_nodes);
    vector<Vertex> vtcsInSPGM = vector<Vertex>(num_nodes);
    vtcsInSPGM[0] = spgmRootV;
    for (int i = 0; i < num_nodes; i++)
    {
        included[i] = false;
    }

    //compute the Conditional Prob Table 
    const CPTfactor& f = ComputeCPT(0, ROOT_VAR, p, verbose);
    outSpgm.SetFactor2Root(spgmRootV, f);

    while (!toProcess.empty())
    {
        Vertex curr = toProcess.front();
        Vertex currInSPGM = vtcsInSPGM[curr];
        toProcess.pop();
        included[curr] = true;

        //find children nodes
        list<Vertex> chds;
        if (verbose)
            cout << "Add children of " << curr << ":\n";
        for (vector < UEdge >::iterator it = mst.begin(); it != mst.end(); it++)
        {
            Vertex v1 = source(*it, fullyConnected);
            Vertex v2 = target(*it, fullyConnected);
            if (v1 == curr || v2 == curr)
            {
                if (verbose)
                    cout << "  process vertex (" << v1 << " - " << v2 << ")\n";
                if (v1 == curr && !included[v2])
                {
                    chds.push_back(v2);
                    if (verbose)
                        cout << "  add " << v2 << endl;
                }
                else if (v2 == curr && !included[v1])
                {
                    chds.push_back(v1);
                    if (verbose)
                        cout << "  add " << v1 << endl;
                }
            }
        }

        Vertex parentInSPGM;
        if (chds.size() <= 1)
        {
            parentInSPGM = currInSPGM;
        }
        else
        {
            //create a sum node, attach it to curr
            parentInSPGM = outSpgm.AddProdNode();
            outSpgm.AddEdge(currInSPGM, parentInSPGM);
        }

        //connect parent to children nodes
        for (list<Vertex>::iterator it = chds.begin(); it != chds.end(); it++)
        {
            //add edges to the nodes which were not already included, add the children to the queue
            Vertex child = *it;
            toProcess.push(child);
            Vertex childInSPGM = outSpgm.AddVNode(child);
            vtcsInSPGM[child] = childInSPGM;
            outSpgm.AddEdge(parentInSPGM, childInSPGM);

            //compute the Conditional Prob Table 
            const CPTfactor& f = ComputeCPT(child, curr, p, verbose);
            outSpgm.SetFactor(childInSPGM, currInSPGM, f);
        }
    }

    outSpgm.m_verbose = verbose;
    outSpgm.Finalize();
    outSpgm.m_verbose = false;
    //    endTime = System.nanoTime();
    //
    //    duration = (endTime - startTime) / 1000000000.0;
    //
    //    if (m_printTimes)
    //    {
    //        printf("computing ChowLiu remaining time: %f sec\n", duration);
    //    }

}

CPTfactor ComputeCPT(int childV, int parentV, const ProbabilitiesStruct& p, bool verbose)
{
    //todo memory leak in creating fators (due to laziness). Not a problem in practice for many cases, but it can become one. 

    int v = childV;
    int w = parentV;
    if (parentV != ROOT_VAR)
    {
        //i j is probability of state i given parent state j 
        Real cpt[2][2];
        const RealMatrix& p00 = p.p00;
        const RealMatrix& p10 = p.p10;
        const RealMatrix& p01 = p.p01;
        const RealMatrix& p11 = p.p11;

        if (p00.GetVal(v, w) == 0 && p10.GetVal(v, w) == 0)
        {
            if (verbose)
                cerr << "too little data to estimate CPT";
            cpt[0][0] = 0.5;
            cpt[1][0] = 0.5;
        }
        else
        {
            Real z = p00.GetVal(v, w) + p10.GetVal(v, w);
            cpt[0][0] = p00.GetVal(v, w) / z;
            cpt[1][0] = p10.GetVal(v, w) / z;
        }

        if (p01.GetVal(v, w) == 0 && p11.GetVal(v, w) == 0)
        {
            if (verbose)
                cerr << "too little data to estimate CPT";
            cpt[0][1] = 0.5;
            cpt[1][1] = 0.5;
        }
        else
        {
            Real z = p01.GetVal(v, w) + p11.GetVal(v, w);
            cpt[0][1] = p01.GetVal(v, w) / z;
            cpt[1][1] = p11.GetVal(v, w) / z;
        }

        if (verbose)
        {
            printf("CPT(%d | %d):[%f, %f; %f, %f]\n", v, w, cpt[0][0], cpt[0][1], cpt[1][0], cpt[1][1]);
        }

        CPTfactor f(2, 2);
        //note the transpose
        f.SetWeight(0, 0, cpt[0][0]);
        f.SetWeight(0, 1, cpt[1][0]);
        f.SetWeight(1, 0, cpt[0][1]);
        f.SetWeight(1, 1, cpt[1][1]);
        return f;
    }
    else
    {
        Real cpt[2];

        if (p.p0[v] == 0 && p.p1[v] == 0)
        {
            cerr << "too little data to estimate CPT";
            cpt[0] = 0.5;
            cpt[1] = 0.5;
        }
        else
        {
            Real z = p.p0[v] + p.p1[v];
            cpt[0] = p.p0[v] / z;
            cpt[1] = p.p1[v] / z;
        }

        if (verbose)
        {
            printf("CPT(%d):[%f, %f]\n", v, cpt[0], cpt[1]);
        }

        CPTfactor f(1, 2);
        f.SetWeight(0, 0, cpt[0]);
        f.SetWeight(0, 1, cpt[1]);
        return f;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

PathInfo FindPath(SPGM& treeSpgm, varToNodeMap map, int varA, int varB, bool verbose)
{
    NodeInfo& na = map[varA];
    NodeInfo& nb = map[varB];
    Vertex currVertex = EMPTY_VERTEX;
    int startVar;
    int otherVar;
    int loopRootVar = -1;
    Vertex otherVertex;
    if (na.depth > nb.depth)
    {
        startVar = varA;
        otherVar = varB;
        currVertex = na.Vertex;
        otherVertex = nb.Vertex;
    }
    else
    {
        startVar = varB;
        otherVar = varA;
        currVertex = nb.Vertex;
        otherVertex = na.Vertex;
    }

    int nvarsInSPGM = map.size();
    vector<bool> visitedVars(nvarsInSPGM);
    for (int i = 0; i < nvarsInSPGM; i++)
        visitedVars[i] = false;

    PathInfo path;

    //first go from currVertex to the root, stopping if we find otherVar
    bool otherVarReached = false;
    int lastVar = -1;
    while (true)
    {
        //if the vertex is a vnode, then save the variable
        const SPGMnode* currSPGMnode = treeSpgm.EditNode(currVertex);
        if (currSPGMnode->GetType() == SPGMnode::vnode)
        {
            Vnode* vn = (Vnode*) currSPGMnode;
            int currVar = vn->GetVar();
            visitedVars[currVar] = true;

            if (lastVar != -1)
            {
                path.varChdMap[currVar] = list<int>();
                path.varChdMap[currVar].push_back(lastVar);
            }
            lastVar = currVar;

            if (currVar == otherVar)
            {
                otherVarReached = true;
                loopRootVar = otherVar;
                break; //top of the path
            }
        }

        //move to the vertex above
        const list<Vertex>& parents = treeSpgm.GetParentVertices(currVertex);
        if (parents.empty())
            break; //we arrived at the root;
        else
        {
            assert(parents.size() == 1);
            Vertex parent = parents.front();

            path.vertChildMap[parent] = list<Vertex>();
            path.vertChildMap[parent].push_back(currVertex);

            currVertex = parent;
            //            vertices.push_front(parent);
        }
    }

    //if we did not reach otherVar, then start from otherVar and go up until a visited var is found
    if (!otherVarReached)
    {
        currVertex = otherVertex;
        lastVar = -1;
        while (true)
        {
            const SPGMnode* currSPGMnode = treeSpgm.EditNode(currVertex);
            if (currSPGMnode->GetType() == SPGMnode::vnode)
            {
                Vnode* vn = (Vnode*) currSPGMnode;
                int currVar = vn->GetVar();

                if (lastVar != -1)
                {
                    //add to map 
                    VarToChildrenVarMap::iterator found = path.varChdMap.find(currVar);
                    if (found == path.varChdMap.end())
                    {
                        list<int> temp = list<int>();
                        temp.push_back(lastVar);
                        path.varChdMap[currVar] = temp;
                    }
                    else
                    {
                        found->second.push_back(lastVar);
                        assert(found->second.size() <= 2);
                    }
                }
                lastVar = currVar;

                if (visitedVars[currVar])
                {
                    loopRootVar = currVar;
                    break;
                }
                else
                {
                    visitedVars[currVar] = true;
                    assert(currVar != startVar);
                }
            }

            //move to the vertex above
            const list<Vertex>& parents = treeSpgm.GetParentVertices(currVertex);
            assert(parents.size() == 1);
            Vertex parent = parents.front();

            VertexToChildrenVerticesMap::iterator found = path.vertChildMap.find(parent);
            if (found == path.vertChildMap.end())
            {
                list<Vertex> temp = list<Vertex>();
                temp.push_back(currVertex);
                path.vertChildMap[parent] = temp;
            }
            else
            {
                list<Vertex>& currAddedChds = found->second;
                if (currAddedChds.empty() || currAddedChds.front() != currVertex)
                {
                    currAddedChds.push_back(currVertex);
                    assert(currAddedChds.size() <= 2);
                }
            }

            currVertex = parent;
        }

        //now go back from currVertex to the root and remove the components from the adj matrix 
        const list<Vertex>& parents = treeSpgm.GetParentVertices(currVertex);
        if (!parents.empty())
        {
            assert(parents.size() == 1);
            Vertex parent = parents.front();
            currVertex = parent;

            while (true)
            {
                //if the vertex is a vnode, then save the variable
                path.vertChildMap.erase(currVertex);
                const SPGMnode* currSPGMnode = treeSpgm.EditNode(currVertex);
                if (currSPGMnode->GetType() == SPGMnode::vnode)
                {
                    Vnode* vn = (Vnode*) currSPGMnode;
                    int currVar = vn->GetVar();
                    visitedVars[currVar] = false;
                    path.varChdMap.erase(currVar);
                }

                //move to the vertex above
                const list<Vertex>& parents = treeSpgm.GetParentVertices(currVertex);
                if (parents.empty())
                    break; //we arrived at the root;
                else
                {
                    assert(parents.size() == 1);
                    Vertex parent = parents.front();
                    currVertex = parent;
                }
            }
        }
    }

    path.rootVar = loopRootVar;

    //compute vector of all traversed vars
    int nvars = 0;
    for (int i = 0; i < nvarsInSPGM; i++)
        if (visitedVars[i])
            nvars++;
    path.m_vars = vector<int>(nvars);
    int ind = 0;
    for (int i = 0; i < nvarsInSPGM; i++)
        if (visitedVars[i])
            path.m_vars[ind++] = i;

    if (verbose)
    {
        cout << "created path between " << varA << " and " << varB << ":\n" << path.ToString();
    }

    return path;

}

void FindMinEdgeInPath(PathInfo& path, Matrix_rowMajor<Real>& mi, bool verbose)
{
    path.minEdge.w = std::numeric_limits<Real>::max();

    for (VarToChildrenVarMap::iterator it = path.varChdMap.begin(); it != path.varChdMap.end(); it++)
    {
        int currVar = it->first;
        const list<int>& vl = it->second;
        for (list<int>::const_iterator it2 = vl.begin(); it2 != vl.end(); it2++)
        {
            int nextVar = *it2;
            Real edgeCost = mi.GetVal(currVar, nextVar);
            if (edgeCost < path.minEdge.w)
            {
                path.minEdge.a = currVar;
                path.minEdge.b = nextVar;
                path.minEdge.w = edgeCost;
            }
        }
    }
    if (verbose)
        printf("min edge in path is (%d,%d) with MI %f\n", path.minEdge.a, path.minEdge.b, path.minEdge.w);
}

VarToChildrenMap ComputeBorderNodes(SPGM& spgm, VertexToChildrenVerticesMap& v2chInTreeSPGM, PathInfo& path, varToNodeMap& treeSpgmVarToNodeMap, bool verbose)
{
    VarToChildrenMap borderNodes;
    //for each variable in the path, computes a list of the vertices that should be its children and are outside of the path. 
    for (int i = 0; i < path.m_vars.size(); i++)
    {
        int currVar = path.m_vars[i];
        borderNodes[currVar] = list<Vertex>();

        NodeInfo& currNodeInfo = treeSpgmVarToNodeMap[currVar];
        list<Vertex>& allChildrenVertices = v2chInTreeSPGM[currNodeInfo.Vertex]; //children also outside path

        //it's a vnode, so it has at most 1 vertex
        if (allChildrenVertices.empty())
            continue;
        assert(allChildrenVertices.size() == 1);

        //take the child of the current vNode who lies in the path
        Vertex childInPath = EMPTY_VERTEX;
        list<Vertex>& childrenInPath = path.vertChildMap[currNodeInfo.Vertex]; //take the vertex descending from the original vertex in the path
        if (!childrenInPath.empty())
            childInPath = childrenInPath.front();

        //take the child vertex 
        Vertex childVertex = allChildrenVertices.front();
        SPGMnode* childNode = spgm.EditNode(childVertex);

        //if it is a product node, then 
        if (childNode->GetType() == SPGMnode::prod && childInPath != EMPTY_VERTEX)
        {
            assert(childVertex == childInPath); //the child in path must be childVertex
            Vertex prodVertex = childVertex;

            //find the children of the prod node. If it is the loop root, it can have two children in the path, otherwise one. 
            //check which children are in the path and which are not. Add the ones that are not to border nodes. 
            list<Vertex>& grandchildrenInPath = path.vertChildMap[prodVertex]; //children in path
            const list<Vertex>& allGrandchildren = v2chInTreeSPGM[prodVertex]; //children outside path
            assert(!allGrandchildren.empty());
            assert(grandchildrenInPath.size() <= 2);

            list<Vertex>& bn = borderNodes[currVar];
            for (list<Vertex>::const_iterator it = allGrandchildren.begin(); it != allGrandchildren.end(); it++)
            {
                Vertex ch = *it;
                bool isInPath = false;
                for (list<Vertex>::const_iterator inPath = grandchildrenInPath.begin(); inPath != grandchildrenInPath.end(); inPath++)
                {
                    Vertex childInPath = *inPath;
                    if (childInPath != EMPTY_VERTEX && ch == childInPath)
                    {
                        isInPath = true;
                        break;
                    }
                }
                //if vertex not in path, add it
                if (!isInPath)
                    bn.push_back(ch);
            }
        }
            //if it is not a prod node it is a vNode. 
        else
        {
            //if vertex not in path, add it
            if (childVertex != childInPath)
            {
                borderNodes[currVar].push_back(childVertex);
            }
        }



    }
    if (verbose)
    {
        cout << "border nodes: \n" << VarToChildrenMapToString(borderNodes);
    }
    return borderNodes;
}

PathInfo CreateAlternativePath(PathInfo& originalPath, bool verbose)
{
    /*takes a directed path from A to B, removes the min edge, adds another edge, creates new directed path with same root*/

    if (verbose)
        cout << "CreateAlternativePath removing edge (" << originalPath.minEdge.a << "," << originalPath.minEdge.b << ")\n";

    int nvars = originalPath.m_vars.size();

    //create undirected adj matrix by adding direct and inverse connections of old path (opposite arrow direction),
    //but remove the min edge and add the new edge to the undirected adj matrix
    VarToChildrenVarMap undirNewAdj = originalPath.varChdMap;
    for (int i = 0; i < nvars; i++)
        undirNewAdj[originalPath.m_vars[i]] = list<int>();
    for (VarToChildrenVarMap::iterator it = originalPath.varChdMap.begin(); it != originalPath.varChdMap.end(); it++)
    {
        int va = it->first;
        const list<int>& vl = it->second;
        for (list<int>::const_iterator it2 = vl.begin(); it2 != vl.end(); it2++)
        {
            int vb = *it2;
            bool isEdgeToRemove = (va == originalPath.minEdge.a && vb == originalPath.minEdge.b) || (va == originalPath.minEdge.b && vb == originalPath.minEdge.a);
            if (!isEdgeToRemove)
            {
                undirNewAdj[vb].push_back(va); //direct connection
                undirNewAdj[va].push_back(vb); //reverse connection
            }
        }
    }
    undirNewAdj[originalPath.edgeToAdd.a].push_back(originalPath.edgeToAdd.b); //new edge
    undirNewAdj[originalPath.edgeToAdd.b].push_back(originalPath.edgeToAdd.a);
    if (verbose)
        cout << "Undirected adj matrix of new path: \n" << VarToChildrenVarMapToString(undirNewAdj);

    //starting from the root, visit the undirected adj matrix and create a directed version
    PathInfo altPath;
    altPath.rootVar = originalPath.rootVar;
    for (int i = 0; i < nvars; i++)
        altPath.varChdMap[originalPath.m_vars[i]] = list<int>();
    map<int, bool> visited;
    for (int i = 0; i < nvars; i++)
        visited[originalPath.m_vars[i]] = false;
    visited[originalPath.rootVar] = true;
    queue<int> toProcess;
    toProcess.push(originalPath.rootVar);
    while (!toProcess.empty())
    {
        int currVar = toProcess.front();
        toProcess.pop();
        const list<int>& neighbors = undirNewAdj[currVar];

        //add as children all connections that were not visited, mark them visited and add them to the queue. 
        for (list<int>::const_iterator it = neighbors.begin(); it != neighbors.end(); it++)
        {
            int neighVar = *it;
            if (!visited[neighVar])
            {
                visited[neighVar] = true;
                toProcess.push(neighVar);
                altPath.varChdMap[currVar].push_back(neighVar);
            }
        }
    }


    if (verbose)
    {
        printf("created alternative path adding edge (%d,%d) to original path:\n", originalPath.edgeToAdd.a, originalPath.edgeToAdd.b);
        //        cout << originalPath.ToString() << "\nPath2:\n" << altPath.ToString() << endl;
        cout << altPath.ToString() << endl;
    }

    return altPath;
}

void AddAlternativePath(SPGM& spgm, PathInfo& originalPath, PathInfo& newPath, VarToChildrenMap& borderNodes, varToNodeMap& treeSpgmVarToNodeMap, map<Vertex, list<PathCoupleInfo> >& sumWeightsMap, bool verbose)
{
    // first, do not insert the path root, because it's in common. 
    //however, check that the path root in the tree SPGM has a sum node as child. If not, create it.
    //use the sum node below the root as currNode

    const NodeInfo& loopRoot = treeSpgmVarToNodeMap[newPath.rootVar];
    Vnode* loopRootNode = (Vnode*) loopRoot.spgmNode;

    //find the child of the loop root
    list<Vertex> rcv = spgm.GetChildrenVertices(loopRoot.Vertex);
    assert(rcv.size() == 1);
    Vertex rootChildVertex = rcv.front();
    const SPGMnode* rootChildNode = spgm.EditNode(rootChildVertex);

    Real originalPathW = originalPath.minEdge.w; //the min weight (in the original path in the MST) that is replaced in the new path
    Real alternativePathW = originalPath.edgeToAdd.w;
    assert(alternativePathW <= originalPathW); //otherwise the MST would not be maximal

    //todo check that each sum weight is assigned to the right child! todo, use assignment through graph edge weights...
    //todo see if it's better to use whole path weight. 
    NodeInfo loopTop;
    if (rootChildNode->GetType() == SPGMnode::sum)
    {
        loopTop.Vertex = rootChildNode->GetVertex();

        //add info for computing sum node weights
        PathCoupleInfo pc;
        pc.original = originalPath;
        pc.new = newPath;
        sumWeightsMap[rootChildNode->GetVertex()].push_back(pc);
    }
    else
    {
        Vertex oldChildVertex = rootChildNode->GetVertex();
        spgm.RemoveEdge(loopRoot.Vertex, oldChildVertex);
        Vertex rootSum = spgm.AddSumNode(); //create new sum node
        spgm.AddEdge(loopRoot.Vertex, rootSum);
        spgm.AddEdge(rootSum, rootChildNode->GetVertex());
        //        spgm.AddSumEdge(rootSum, rootChildNode->GetVertex(), originalPathW);
        loopTop.Vertex = rootSum;

        // add info for computing sum node weights
        list<PathCoupleInfo> pciList;

        PathCoupleInfo pcOrig;
        pcOrig.original = originalPath;
        pcOrig.new = originalPath; //no alternative path here 
        pciList.push_back(pcOrig);

        PathCoupleInfo pc;
        pc.original = originalPath;
        pc.new = newPath;
        pciList.push_back(pc);

        sumWeightsMap[rootSum] = pciList;
    }

    //Attach the nodes below to loopTop
    list<int>& rootChdsInside = newPath.varChdMap[loopRoot.var]; // children of root inside loop;
    list<Vertex>& rootChdsOutside = borderNodes[loopRoot.var]; //children Of  Root outside loop

    queue<NodeInfo> toProcess;

    if ((rootChdsInside.size() + rootChdsOutside.size()) > 1)
    {
        //if n children > 1, create a product
        Vertex prod = spgm.AddProdNode();
        //        spgm.AddSumEdge(loopTop.Vertex, prod, alternativePathW);
        spgm.AddEdge(loopTop.Vertex, prod);

        loopTop.Vertex = prod;
        //connect the vertex to all the children outside 
        for (list<Vertex>::iterator it = rootChdsOutside.begin(); it != rootChdsOutside.end(); it++)
            spgm.AddEdge(loopTop.Vertex, *it);
        //create new children inside and connect to them, then add them to process
        for (list<int>::iterator it = rootChdsInside.begin(); it != rootChdsInside.end(); it++)
        {
            NodeInfo newNode;
            newNode.var = *it;
            newNode.Vertex = spgm.AddVNode(newNode.var);
            spgm.AddEdge(loopTop.Vertex, newNode.Vertex);
            toProcess.push(newNode);
        }
    }
    else
    {
        //there must be only 1 child, inside the loop 
        assert(rootChdsInside.size() == 1 && rootChdsOutside.empty());
        int childVar = rootChdsInside.front();
        NodeInfo newNode;
        newNode.var = childVar;
        newNode.Vertex = spgm.AddVNode(newNode.var);
        //        spgm.AddSumEdge(loopTop.Vertex, newNode.Vertex, alternativePathW);
        spgm.AddEdge(loopTop.Vertex, newNode.Vertex);
        toProcess.push(newNode);
    }

    while (!toProcess.empty())
    {
        //take the first vertex in the queue
        NodeInfo currNode = toProcess.front();
        toProcess.pop();
        int currVar = currNode.var;

        const list<int>& chdsInside = newPath.varChdMap[currVar]; // children of currNode inside loop;
        const list<Vertex>& chdsOutside = borderNodes[currVar]; //children Of currNode outside loop

        if (chdsInside.size() + chdsOutside.size() > 1)
        {
            //if n children > 1, create a product
            Vertex prod = spgm.AddProdNode();
            spgm.AddEdge(currNode.Vertex, prod);
            currNode.Vertex = prod;
        }

        //connect the vertex to all the children outside 
        for (list<Vertex>::const_iterator it = chdsOutside.begin(); it != chdsOutside.end(); it++)
            spgm.AddEdge(currNode.Vertex, *it);
        //create new children inside and connect to them, then add them to process
        for (list<int>::const_iterator it = chdsInside.begin(); it != chdsInside.end(); it++)
        {
            NodeInfo newNode;
            newNode.var = *it;
            newNode.Vertex = spgm.AddVNode(newNode.var);
            spgm.AddEdge(currNode.Vertex, newNode.Vertex);
            toProcess.push(newNode);
        }
    }
}

void AddFactors(SPGM& spgm, ProbabilitiesStruct& probs, bool verbose)
{
    //adds factors from each node to each parent, based on probability table probs. 
    //Note that the same factor is used for same parent variables, therefore the factors are shared whenever they connect the same variables. 
    typedef std::pair<int, int> varEdge;
    typedef map< varEdge, CPTfactor> edgeFactorMap;
    edgeFactorMap efm;
    spgm.ComputeTopOrder();
    spgm.ComputeVParents();

    //collect information of all edges in the graph
    int N = spgm.GetNNodes();
    for (int i = 0; i < N; i++)
    {
        const SPGMnode* n = spgm.GetNodeInOrder(i);
        if (n->GetType() != SPGMnode::vnode)
            continue;
        const Vnode* vn = (const Vnode*) n;

        int currVar = vn->GetVar();
        const std::set<int>& vpas = spgm.GetVParents(vn->GetVertex());
        if (vpas.empty())
        {
            varEdge edge;
            edge.first = ROOT_VAR;
            edge.second = currVar;
            efm[edge] = CPTfactor(); //for now just collect which edges are to add, then we add the factor later
        }
        else
        {
            for (set<int>::iterator it = vpas.begin(); it != vpas.end(); it++)
            {
                varEdge edge;
                edge.first = *it;
                edge.second = currVar;
                efm[edge] = CPTfactor(); //for now just collect which edges are to add, then we add the factor later
            }
        }
    }

    if (verbose)
    {
        cout << "Adding conditional factors for the following variable pairs:\n";
        for (edgeFactorMap::iterator it = efm.begin(); it != efm.end(); it++)
            cout << "    var " << it->first.second << " given " << it->first.first << ",\n";
        cout << "\n";
    }

    //adding the factors
    for (edgeFactorMap::iterator it = efm.begin(); it != efm.end(); it++)
    {
        int childVar = it->first.second;
        int parentVar = it->first.first;
        const CPTfactor& f = ComputeCPT(childVar, parentVar, probs, verbose);
        it->second = f;
    }

    //settting the factors to each node
    for (int i = 0; i < N; i++)
    {
        const SPGMnode* n = spgm.GetNodeInOrder(i);
        Vertex currVertex = n->GetVertex();
        if (n->GetType() != SPGMnode::vnode)
            continue;
        const Vnode* vn = (const Vnode*) n;

        int currVar = vn->GetVar();
        //        const std::set<int>& vpas = spgm.GetVParents(vn->GetVertex());
        const std::set<Vertex>& vpaVertices = spgm.GetVParentVertices(vn->GetVertex());
        for (set<Vertex>::iterator it = vpaVertices.begin(); it != vpaVertices.end(); it++)
        {
            Vertex paVert = *it;
            if (paVert == ROOT_VERTEX)
            {
                varEdge edge;
                edge.first = ROOT_VAR;
                edge.second = currVar;
                const CPTfactor& f = efm[edge];
                spgm.SetFactor2Root(currVertex, f);
            }
            else
            {
                varEdge edge;
                Vnode* paVn = (Vnode*) spgm.EditNode(*it);
                edge.first = paVn->GetVar();
                edge.second = currVar;
                const CPTfactor& f = efm[edge];
                spgm.SetFactor(currVertex, *it, f);
            }
        }
    }
}

void ApplySumWeights(SPGM& spgm, map<Vertex, list<PathCoupleInfo> >& wmap)
{
//COMPUTE WEIGHTS FROM PATH INFOS 
//weight is proportional to the sum of log probabilities of the best tree at child c
// therefore, only the trees in the loop from c to the lower node are considered
 
    for (map<Vertex, list<PathCoupleInfo> >::iterator itw = wmap.begin(); itw != wmap.end(); itw++)
    {
        int vertex = itw->first;
        SumNode* sn = (SumNode*) spgm.EditNode(vertex);
        const list<PathCoupleInfo>& paths = itw->second;

        //FIRST, LOOK FOR LONGEST ORIGINAL PATH - call it T
        PathInfo maxp;
        for (list<PathCoupleInfo>::const_iterator itp = paths.begin(); itp!= paths.end(); itp++)
        {
            PathInfo& path = *itp;
            if (path.m_vars.size() > maxp.m_vars.size())
            {
                maxp = path;
            }
        }

        //THEN, FOR EVERY OTHER TREE, COMPUTE TOTAL WEIGHT AS WEIGHT IN LOOP to B + WEIGHT FROM B to BASE OF T 
        //FACTORS AND LOG P OF FACTORS SHOULD HAVE BEEN COMPUTED
        vector<Real> wvec(paths.size());
        Real sum = 0;
        int i=0;
        for (list<PathCoupleInfo>::const_iterator itp = paths.begin(); itp!= paths.end(); itp++)
        {
            PathInfo& pathCouple = *itp;
            Real pathLL = pathCouple.new.LL + maxp.LL - pathCouple.original.LL;
            sum += pathLL;
            wvec[i++] = pathLL;
        }

        //NORMALIZE
        for (int j=0; j<wvec.size(); j++)
            wvec[j++] = *it / sum;

        sn->SetWeights(wvec);
    }
}

void LearnMixMaxSpanningTreesSPGM(SPGM& outSpgm, const DataMatrix& data, const Params& p)
{
    std::vector<Real> emptyVec = std::vector<Real>();
    LearnMixMaxSpanningTreesSPGM(outSpgm, data, emptyVec, p);
}

void LearnMixMaxSpanningTreesSPGM(SPGM& outSpgm, const DataMatrix& data, std::vector<Real> weights, const Params& p)
{
    //• [S, unusedEdges] <-CreateChowLiuSPG(data) // S contains PRODUCT NODES for all splits. 
    int NaddEdges = p.nInsertions;
    Real diricheletPrior = p.diricheletPrior;
    Real lambda = p.lambda;
    bool verbose = p.verbose;
    if (verbose)
    {
        printf("###################################################\nLearnMixMaxSpanningTreesSPGM with %d NaddEdges and diricheletPrior %f\n", NaddEdges, diricheletPrior);
    }




    Matrix_rowMajor<Real> mi;
    ProbabilitiesStruct probs;

    ComputeMutualInfo_binary(data, weights, mi, probs, diricheletPrior, false);

    if (verbose)
    {
        cout << "Mutual Information is:\n" << mi.ToString();
    }

    outSpgm = SPGM();
    ChowLiuTree(outSpgm, mi, probs, false);
    if (verbose)
    {
        printf("learned ChowLiuTree (and saved .dot file in folder LearnMixMaxSpanningTreesSPGM_out)\n");
        outSpgm.WriteGraphViz("LearnMixMaxSpanningTreesSPGM_out/ChowLiuTree");
    }

    //fill information in varToNodeMapInTreeSPGM and vertexToChdInTreeSPGM
    varToNodeMap varToNodeMapInTreeSPGM;
    VertexToChildrenVerticesMap vertexToChdInTreeSPGM;
    vector<const SPGMnode*> invTopOrderTree = outSpgm.GetInvTopOrder();
    for (int i = 0; i < invTopOrderTree.size(); i++)
    {
        const SPGMnode* np = invTopOrderTree[i];
        NodeInfo n;
        n.depth = invTopOrderTree.size() - i - 1;
        n.Vertex = np->GetVertex();
        n.var = -1;
        n.spgmNode = np;
        vertexToChdInTreeSPGM[n.Vertex] = outSpgm.GetChildrenVertices(n.Vertex);

        if (np->GetType() == SPGMnode::vnode)
        {
            Vnode* vnp = (Vnode*) np;
            n.var = vnp->GetVar();
        }
        varToNodeMapInTreeSPGM[n.var] = n;
    }

    //find list of edges not included in MST and corresponding paths in MST connecting the two elements of the edge
    int num_nodes = mi.nCols();
    assert(mi.nCols() == mi.nRows());
    list<PathInfo> pathInfosToAdd;
    if (verbose)
        printf("\n######################################################\nProcessing potential edges to add:\n\n");
    if (NaddEdges > 0)
    {
        for (int i = 0; i < num_nodes; i++)
        {
            for (int j = i + 1; j < num_nodes; j++)
            {
                PathInfo path = FindPath(outSpgm, varToNodeMapInTreeSPGM, i, j, verbose);
                assert(path.m_vars.size() > 1);
                bool isInMst = path.m_vars.size() == 2;
                if (!isInMst)
                {
                    WeightedEdgeInfo wei;
                    wei.a = i;
                    wei.b = j;
                    wei.w = mi.GetVal(i, j);
                    path.edgeToAdd = wei;
                    pathInfosToAdd.push_back(path);
                    if (verbose)
                        printf("Include this path\n\n");
                }
                else if (verbose)
                    printf("Edge (%d,%d) already in spgm\n\n", i, j);
            }
        }
    }

    /* solve a knapsack problem for choosing which paths to add. 
       the value of each path is the weight of the added edge
       the cost of each path is the path length */
    list<PathInfo> selectedPathInfos;
    if (lambda < 0)
    {
        cerr << "wrong lambda (must be >=0 ): " << lambda << endl;
        exit(1);
    }
    if (isinf(lambda))
    {
        //just take NaddEdges paths in order of decreasing weights 
        pathInfosToAdd.sort();
        pathInfosToAdd.reverse(); //sort in descending order
        if (verbose)
        {
            cout << "ordered paths weights: ";
            for (list<PathInfo>::iterator it = pathInfosToAdd.begin(); it != pathInfosToAdd.end(); it++)
            {
                cout << (*it).edgeToAdd.w << ", ";
            }
            cout << endl;
        }

        int nAdded = 0;
        for (list<PathInfo>::iterator it = pathInfosToAdd.begin(); it != pathInfosToAdd.end(); it++)
        {
            if (nAdded >= NaddEdges)
                break;
            selectedPathInfos.push_back(*it);
            int edgesToAdd = (*it).m_vars.size() - 1;
            assert(edgesToAdd > 0);
            nAdded += edgesToAdd;
        }
    }
    else
    {
        vector<Real> values(pathInfosToAdd.size());
        vector<int> costs(pathInfosToAdd.size());
        int i = 0;
        for (list<PathInfo>::iterator it = pathInfosToAdd.begin(); it != pathInfosToAdd.end(); it++)
        {
            values[i] = pow((*it).edgeToAdd.w, lambda); //the value is defined as the edge weight that we add
            int theCost = ((*it).m_vars.size() - 1); //the cost is the number of edges that we create 
            assert(theCost > 0);
            costs[i] = theCost;
            i++;
        }
        vector<bool> chosen;
        KnapSack(p.nInsertions, costs, values, chosen);

        i = 0;
        for (list<PathInfo>::iterator it = pathInfosToAdd.begin(); it != pathInfosToAdd.end(); it++)
        {
            if (chosen[i])
                selectedPathInfos.push_back(*it);
            i++;
        }
    }




    /* add all the selected path infos */
    map<Vertex, list<PathCoupleInfo> > sumWeightsMap;
    if (verbose)
        printf("\n######################################################\nAdd all the selected path infos:\n\n");
    int i = 0;
    for (list<PathInfo>::iterator it = selectedPathInfos.begin(); it != selectedPathInfos.end(); it++)
    {
        PathInfo& originalPath = *it;
        if (verbose)
            cout << "ADDING EDGE (" << originalPath.edgeToAdd.a << "," << originalPath.edgeToAdd.b << "), weight " << originalPath.edgeToAdd.w << " TO PATH:\n" << originalPath.ToString();

        FindMinEdgeInPath(originalPath, mi, verbose);
        VarToChildrenMap borderNodes = ComputeBorderNodes(outSpgm, vertexToChdInTreeSPGM, originalPath, varToNodeMapInTreeSPGM, verbose);
        PathInfo newPath = CreateAlternativePath(originalPath, verbose);
        AddAlternativePath(outSpgm, originalPath, newPath, borderNodes, varToNodeMapInTreeSPGM, sumWeightsMap, verbose);

        if (verbose)
        {
            stringstream filename;
            filename << "LearnMixMaxSpanningTreesSPGM_out/LearnMixMaxSpanningTreesSPGM_insertion" << i;
            printf("Writing SPGM at insertion %d on file:\n %s.dot\n\n", i, filename.str().c_str());
            outSpgm.WriteGraphViz(filename.str());
        }
        i++;
    }

    outSpgm.m_verbose = verbose;
    outSpgm.ClearMessages();
    ApplySumWeights(outSpgm, sumWeightsMap);
    AddFactors(outSpgm, probs, verbose);
    outSpgm.Finalize();
    outSpgm.ComputeScope();
    outSpgm.m_verbose = false;
    //    fprintf("LearnMixMaxSpanningTreesSPGM done. \nVnode edges in SPGM %d, in mixture of trees %d.\nSaved %d edges (%.1f percent). \n",usedEdges,totalEdges,totalEdges-usedEdges, 100.f*((Real)usedEdges/(Real)totalEdges));
    if (verbose)
        printf("LearnMixMaxSpanningTreesSPGM done. Added %d paths. \n", (int) selectedPathInfos.size());

}

void LearnMixEM(const Params& p, Real* out_trainLL , Real* out_validLL, Real* out_testLL)
{
    printf("\nLearnMixEm()\n");
    //    printf("\n############\nLearnMixEm() with params:\n%s\n############\n", p.ToString().c_str());

    string trainFile = p.dataFolder + p.filename + ".ts.data";
    //        string trainFile = "testData/test.ts";
    DataMatrix train;
    train.ReadFromCSV(trainFile);
    printf("loaded %s train with %d vars %d samples \n", p.filename.c_str(), (int) train.nVars(), (int) train.nData());

    string validFile = p.dataFolder + p.filename + ".valid.data";
    DataMatrix valid;
    valid.ReadFromCSV(validFile);
    printf("loaded valid data with %d vars %d samples \n", (int) valid.nVars(), (int) valid.nData());
    
    string testFile = p.dataFolder + p.filename + ".test.data";
    //        string testFile = "testData/test.ts";
    DataMatrix test;
    test.ReadFromCSV(testFile);
    printf("loaded test data with %d vars %d samples \n", (int) test.nVars(), (int) test.nData());

    SPGM_mix spgmMix(p);

    if (p.init == Params::init_kmeans)
        spgmMix.Init_Kmeans(train, p);
    else if (p.init == Params::init_rand)
        spgmMix.Init_Random(train, p);

    printf("initial spgmMix has %d messages\n", spgmMix.GetNMessages());

    Real testLL = spgmMix.MeanLL(test);
    printf("\nINITAL test LL: %.12f\n", testLL);

    spgmMix.m_verbose = false;
    spgmMix.EM_struct(train, test, p);
    printf("spgmMix has %d messages\n", spgmMix.GetNMessages());
    spgmMix.EM_params(train, test, p);
    printf("spgmMix has %d messages\n", spgmMix.GetNMessages());

    //    spgmMix.EM_struct(train, test, p);
    //    printf("spgmMix has %d messages\n",spgmMix.GetNMessages());
    //    spgmMix.EM_params(train, test, p);
    //    printf("spgmMix has %d messages\n",spgmMix.GetNMessages());

    testLL = spgmMix.MeanLL(test);
    printf("\nFINAL test LL %.12f \n", testLL);

    if (out_testLL != NULL)
        *out_testLL = testLL;
    if (out_trainLL != NULL)
        *out_trainLL = spgmMix.MeanLL(train);
    if (out_validLL != NULL)
        *out_validLL = spgmMix.MeanLL(valid);
}