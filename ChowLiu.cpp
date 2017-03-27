
#include "LearnSPGM.h"
#include <boost/graph/kruskal_min_spanning_tree.hpp>

using namespace std;

using namespace boost;
typedef adjacency_list < vecS, vecS, undirectedS, no_property, property < edge_weight_t, Real > > Ugraph;
typedef graph_traits < Ugraph >::edge_descriptor UEdge;
typedef std::pair<int, int> E;

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

        const vector<byte>& dataRow = data.ByRow().GetRow(i);

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
        p0[v] = c0[v] / Nsamples;
        p1[v] = c1[v] / Nsamples;

        if (diricheletPrior > 0)
        {
            Real z = p0[v] + p1[v] + diricheletPrior;
            p0[v] = (p0[v] + diricheletPrior/2) / z;
            p1[v] = (p1[v] + diricheletPrior/2) / z;
        }

        h[v] = -p0[v] * SafeLog(p0[v]) - p1[v] * SafeLog(p1[v]);

        for (int w = v + 1; w < Nnodes; w++)
        {
            Real z = c00.GetVal(v, w) + c01.GetVal(v, w) + c10.GetVal(v, w) + c11.GetVal(v, w);
            p00.SetVal(v, w, (c00.GetVal(v, w)) / z);
            p01.SetVal(v, w, (c01.GetVal(v, w)) / z);
            p10.SetVal(v, w, (c10.GetVal(v, w)) / z);
            p11.SetVal(v, w, (c11.GetVal(v, w)) / z);

            if (diricheletPrior > 0)
            {
                Real z = p00.GetVal(v, w) + p01.GetVal(v, w) + p10.GetVal(v, w) + p11.GetVal(v, w) + diricheletPrior;
                p00.SetVal(v, w, (p00.GetVal(v, w) + diricheletPrior/4) / z);
                p01.SetVal(v, w, (p01.GetVal(v, w) + diricheletPrior/4) / z);
                p10.SetVal(v, w, (p10.GetVal(v, w) + diricheletPrior/4) / z);
                p11.SetVal(v, w, (p11.GetVal(v, w) + diricheletPrior/4) / z);
            }

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
            Real mi = h[v] + h[w] - h_vw;


            if (mi < 0)
            {
                if (verbose)
                    cerr << "ERROR: negative Mutual Information? Can be due to numerical errors\n";
                mi = 0;
            }

            outputMI.SetVal(v, w, mi);

            if (isnan(mi))
            {
                cerr << "there is a NaN";
            }
            if (verbose)
            {
                printf("Mutual information(%d,%d) = %f\n", v, w, mi);
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
    int nInsertedEdges = 0;
    for (int i = 0; i < num_nodes; i++)
    {
        for (int j = i + 1; j < num_nodes; j++)
        {
            Real val = -mi.GetVal(i, j);
            //            assert(val <= 0);
            if (val != 0)
            {
                edge_array[nInsertedEdges] = E(i, j);
                edge_weights[nInsertedEdges] = val;
                //            cout<<"W"<<edge_weights[ind]<<endl;
                nInsertedEdges++;
            }
        }
    }
    assert(nInsertedEdges <= edge_array.size() && nInsertedEdges <= edge_weights.size());
    edge_array.resize(nInsertedEdges);
    edge_weights.resize(nInsertedEdges);

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

    //find the connected components

    //create SPGM representing the mst
    vector<bool> included = vector<bool>(num_nodes);
    vector<Vertex> vtcsInSPGM = vector<Vertex>(num_nodes);
    for (int i = 0; i < num_nodes; i++)
    {
        included[i] = false;
    }
    int lastRootVar = 0;
    list<Vertex> rootVertices;

    while (true)
    {
        //find first non included var and use it as new root
        int rootVarForTree = -1;
        for (int i = lastRootVar; i < num_nodes; i++)
        {
            if (!included[i])
            {
                rootVarForTree = i;
                lastRootVar = rootVarForTree;
                break;
            }
        }
        if (rootVarForTree == -1)
            break; //all variables have been included

        //create subtree starting from non included var
        Vertex spgmRootV = outSpgm.AddVNode(rootVarForTree);
        rootVertices.push_back(spgmRootV);

        if (verbose)
            cout << "\n    #### SUB-ROOT var: " << lastRootVar << " vertex " << spgmRootV << endl;

        queue<Vertex> toProcess;
        toProcess.push(rootVarForTree);

        //visit nodes in breadth first search and at each node add its children in the SPGM 

        vtcsInSPGM[rootVarForTree] = spgmRootV;


        //compute the Conditional Prob Table 
        const CPTfactor& f = ComputeCPT(rootVarForTree, ROOT_VAR, p, verbose);
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
    }

    //join all the subtrees in a product root node (if there are more than one)
    if (verbose)
        cout << "adding root prod node with children edges:";

    if (rootVertices.size() > 1)
    {
        Vertex rootProd = outSpgm.AddProdNode();
        for (list<Vertex>::iterator it = rootVertices.begin(); it != rootVertices.end(); it++)
        {
            outSpgm.AddEdge(rootProd, *it);
            if (verbose)
                cout << *it << ",";
        }
        if (verbose)
            cout << endl;
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
