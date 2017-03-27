
#include "LearnSPGM.h"
#include <boost/graph/kruskal_min_spanning_tree.hpp>

using namespace std;

using namespace boost;
typedef adjacency_list < vecS, vecS, undirectedS, no_property, property < edge_weight_t, Real > > Ugraph;
typedef graph_traits < Ugraph >::edge_descriptor UEdge;
typedef std::pair<int, int> E;

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

string PathInfo::ToString()
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

bool PathInfo::operator<(const PathInfo& right) const
{
    return edgeToAdd.w < right.edgeToAdd.w;
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

void AddAlternativePath(SPGM& spgm, PathInfo& originalPath, PathInfo& newPath, VarToChildrenMap& borderNodes, varToNodeMap& treeSpgmVarToNodeMap, map<Vertex, list<Real> >& sumWeightsMap, bool verbose)
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

        //add a weight for the new path 
        sumWeightsMap[rootChildNode->GetVertex()].push_back(alternativePathW);
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

        // add sum node weights
        list<Real> newW;
        newW.push_back(originalPathW);
        newW.push_back(alternativePathW);
        sumWeightsMap[rootSum] = newW;
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

void ApplySumWeights(SPGM& spgm, map<Vertex, list<Real> >& wmap)
{
    for (map<Vertex, list<Real> >::iterator itw = wmap.begin(); itw != wmap.end(); itw++)
    {
        int vertex = itw->first;
        SumNode* sn = (SumNode*) spgm.EditNode(vertex);
        const list<Real>& wl = itw->second;
        vector<Real> wvec(wl.size());
        Real sum = 0;
        for (list<Real>::const_iterator it = wl.begin(); it != wl.end(); it++)
            sum += *it;
        int j = 0;
        for (list<Real>::const_iterator it = wl.begin(); it != wl.end(); it++)
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
    MixMax_auxStruct auxStruct = FindCandidatePaths(outSpgm, data, weights, p);
    auxStruct.selectedPathInfos = SelectPaths(auxStruct.pathInfosToAdd, p, false);
    AddPaths(outSpgm, auxStruct, p);
}

MixMax_auxStruct FindCandidatePaths(SPGM& outSpgm, const DataMatrix& data, std::vector<Real> weights, const Params& p)
{
    MixMax_auxStruct auxStruct;

    //• [S, unusedEdges] <-CreateChowLiuSPG(data) // S contains PRODUCT NODES for all splits. 
    int NaddEdges = p.nInsertions;
    Real diricheletPrior = p.diricheletPrior;

    bool verbose = p.verbose;
    if (verbose)
    {
        printf("###################################################\nLearnMixMaxSpanningTreesSPGM with %d NaddEdges and diricheletPrior %f\n", NaddEdges, diricheletPrior);
    }

    Matrix_rowMajor<Real>& mi = auxStruct.mi;

    ProbabilitiesStruct& probs = auxStruct.probs;

    ComputeMutualInfo_binary(data, weights, mi, probs, diricheletPrior, false);

    //collect all possible weight additions, ordered by weight in ascending order
    int N = mi.nRows();
    int ind = 0;
    int nTotalDifferentEdges = -N + (N * (N + 1)) / 2;
    vector<WeightedEdgeInfo> allEdgesToProcess(nTotalDifferentEdges);
    for (int r = 0; r < N; r++)
    {
        const vector<Real>& row = mi.GetRow(r);
        for (int c = r + 1; c < N; c++)
        {
            WeightedEdgeInfo wei;
            wei.a = r;
            wei.b = c;
            wei.w = row[c];
            allEdgesToProcess[ind++] = wei;
        }
    }
    assert(ind == nTotalDifferentEdges);
    sort(allEdgesToProcess.begin(), allEdgesToProcess.end()); //sorted in ASCENDING oder

    if (p.takeNbestMi > 0)
    {
        //put to 0 the non K best entries of MI
        //find the element after which 
        int nNonzeroMi = p.takeNbestMi;
        if (nNonzeroMi > allEdgesToProcess.size())
            nNonzeroMi = allEdgesToProcess.size();

        Real thresholdW = allEdgesToProcess[allEdgesToProcess.size() - nNonzeroMi].w;
        
        if (verbose)
            cout<<"putting to 0 edges below "<<thresholdW<<endl;

        //put to zero mi below thresholdW 
        for (int r = 0; r < N; r++)
        {
            for (int c = r + 1; c < N; c++)
            {
                Real& w = mi.EditVal(r, c);
                if (w < thresholdW)
                {
                    w = 0;
                    mi.SetVal(c, r, 0);
                }
            }
        }

        //also in allEdgesToProcess
        for (int i = 0; i < allEdgesToProcess.size(); i++)
            if (allEdgesToProcess[i].w < thresholdW)
                allEdgesToProcess[i].w = 0;
    }

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
    varToNodeMap& varToNodeMapInTreeSPGM = auxStruct.varToNodeMapInTreeSPGM;
    VertexToChildrenVerticesMap& vertexToChdInTreeSPGM = auxStruct.vertexToChdInTreeSPGM;

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

    //cut the vector of edges to add to maximum size
    int nToAdd = p.maxNAddedPaths;
    if (nToAdd > allEdgesToProcess.size())
        nToAdd = allEdgesToProcess.size();
    vector<WeightedEdgeInfo> edgesToProcess(allEdgesToProcess.end() - nToAdd, allEdgesToProcess.end());
    assert(edgesToProcess.size() == nToAdd);

    //find list of edges not included in MST and corresponding paths in MST connecting the two elements of the edge
    assert(mi.nCols() == mi.nRows());
    list<PathInfo>& pathInfosToAdd = auxStruct.pathInfosToAdd;
    if (verbose)
        printf("\n######################################################\nProcessing potential edges to add:\n\n");

    for (int i = 0; i < edgesToProcess.size(); i++)
    {
        const WeightedEdgeInfo wei = edgesToProcess[i];
        if(wei.w < 0)
            cerr<<"ERROR: negative mutual information\n";
        if (wei.w <= 0)
        {
            continue;
            //            cerr << "edge with weight 0 should not be included!";
            //            assert(false);
            //            exit(1);
        }
        if (verbose)
        {
            cout << "PROCESSING EDGE (" << wei.a << "," << wei.b << "), weight " << wei.w << "\n";
        }

        PathInfo path = FindPath(outSpgm, varToNodeMapInTreeSPGM, wei.a, wei.b, verbose);
        assert(path.m_vars.size() > 1);
        bool isInMst = path.m_vars.size() == 2;
        if (!isInMst)
        {
            path.edgeToAdd = wei;
            pathInfosToAdd.push_back(path);
            if (verbose)
                printf("Include this path\n\n");
        }
        else if (verbose)
            printf("Edge (%d,%d) already in spgm\n\n", wei.a, wei.b);

    }

    return auxStruct;
}

list<PathInfo> SelectPaths(list<PathInfo>& pathInfosToAdd, const Params& p, bool multiplyPerNmix)
{
    /* solve a knapsack problem for choosing which paths to add. 
       the value of each path is the weight of the added edge
       the cost of each path is the path length */
    bool verbose = p.verbose;
    Real lambda = p.lambda;
    int NaddEdges = p.nInsertions;
    if (multiplyPerNmix)
        NaddEdges *= p.K;

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
            values[i] = pow((*it).edgeToAdd.w + 1, lambda); //the value is defined as the edge weight that we add
            int theCost = ((*it).m_vars.size() - 1); //the cost is the number of edges that we create 
            assert(theCost > 0);
            costs[i] = theCost;
            i++;
        }
        vector<bool> chosen;
        KnapSack(NaddEdges, costs, values, chosen);

        i = 0;
        for (list<PathInfo>::iterator it = pathInfosToAdd.begin(); it != pathInfosToAdd.end(); it++)
        {
            if (chosen[i])
                selectedPathInfos.push_back(*it);
            i++;
        }
    }

    return selectedPathInfos;
}

void AddPaths(SPGM& outSpgm, MixMax_auxStruct& auxStruct, const Params& p)
{
    bool verbose = p.verbose;

    /* add all the selected path infos */
    map<Vertex, list<Real> > sumWeightsMap;
    list<PathInfo>& selectedPathInfos = auxStruct.selectedPathInfos;
    VertexToChildrenVerticesMap& vertexToChdInTreeSPGM = auxStruct.vertexToChdInTreeSPGM;
    varToNodeMap& varToNodeMapInTreeSPGM = auxStruct.varToNodeMapInTreeSPGM;

    if (verbose)
        printf("\n######################################################\nAdd all the selected path infos:\n\n");
    int i = 0;
    for (list<PathInfo>::iterator it = selectedPathInfos.begin(); it != selectedPathInfos.end(); it++)
    {
        PathInfo& originalPath = *it;
        if (verbose)
            cout << "ADDING EDGE (" << originalPath.edgeToAdd.a << "," << originalPath.edgeToAdd.b << "), weight " << originalPath.edgeToAdd.w << " TO PATH:\n" << originalPath.ToString();

        FindMinEdgeInPath(originalPath, auxStruct.mi, verbose);
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
    AddFactors(outSpgm, auxStruct.probs, verbose);
    outSpgm.Finalize();
    outSpgm.ComputeScope();
    outSpgm.m_verbose = false;
    //    fprintf("LearnMixMaxSpanningTreesSPGM done. \nVnode edges in SPGM %d, in mixture of trees %d.\nSaved %d edges (%.1f percent). \n",usedEdges,totalEdges,totalEdges-usedEdges, 100.f*((Real)usedEdges/(Real)totalEdges));
    if (verbose)
        printf("LearnMixMaxSpanningTreesSPGM done. Added %d paths. \n", (int) selectedPathInfos.size());

}

void LearnMixEM(const Params& p, Real* out_trainLL, Real* out_validLL, Real* out_testLL)
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

    Real validLL = spgmMix.MeanLL(valid);
    if(p.plotTestLL)
        printf("\nINITIAL valid LL %.5f (testLL %.5f)\n", validLL, spgmMix.MeanLL(test));
    else
        printf("\nINITIAL valid LL %.5f \n", validLL);

    spgmMix.m_verbose = false;
    spgmMix.EM_struct(train, valid, p, &test);
    printf("spgmMix has %d messages\n", spgmMix.GetNMessages());
    spgmMix.EM_params(train, valid, p, &test);
    printf("spgmMix has %d messages\n", spgmMix.GetNMessages());

    //    spgmMix.EM_struct(train, test, p);
    //    printf("spgmMix has %d messages\n",spgmMix.GetNMessages());
    //    spgmMix.EM_params(train, test, p);
    //    printf("spgmMix has %d messages\n",spgmMix.GetNMessages());

    Real testLL = spgmMix.MeanLL(test);
    printf("\nFINAL test LL %.12f \n", testLL);

    if (out_testLL != NULL)
        *out_testLL = testLL;
    if (out_trainLL != NULL)
        *out_trainLL = spgmMix.MeanLL(train);
    if (out_validLL != NULL)
        *out_validLL = spgmMix.MeanLL(valid);
}
