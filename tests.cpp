

#include "tests.h"
#include "LearnSPGM.h"
#include <time.h>



using namespace std;


#include "Eigen/Core"
#include <Eigen/Dense>

void TestEigen(int nSamples, int nVals)
{
    //    using namespace Eigen;
    //      ArrayXXf  m(2,2);
    //  
    //  // assign some values coefficient by coefficient
    //  m(0,0) = 1.0; m(0,1) = 2.0;
    //  m(1,0) = 3.0; m(1,1) = m(0,1) + m(1,0);
    //  
    //  // print values to standard output
    //  cout << m << endl << endl;
    // 
    //  // using the comma-initializer is also allowed
    //  m << 1.0,2.0,
    //       3.0,4.0;
    //     
    //  // print values to standard output
    //  cout << m << endl;

    using namespace Eigen;

    MatrixXd data(nVals, nSamples);

    Array2Xd in(2, nSamples);
    MatrixX2d w(2, 2);
    Array2Xd out(2, nSamples);

    for (int i = 0; i < nSamples; i++)
    {
        in(0, i) = SPN_Rand();
        in(1, i) = SPN_Rand();
    }

    w(0, 0) = 1;
    w(1, 0) = 0;
    w(0, 1) = 0;
    w(1, 1) = 1;

    if (nSamples < 20)
    {
        cout << "w \n" << w << endl;
        cout << "in \n" << in << endl;
    }

    clock_t t1, t2;

    t1 = clock();

    Matrix2Xd reducedData(2, nSamples);
    reducedData.row(0) = data.row(0);
    reducedData.row(1) = data.row(1);

    //    minRow = (in.row(0)).min((in.row(1)));
    ArrayXd maxRow = in.colwise().minCoeff();

    if (nSamples < 20)
    {
        cout << "maxRow: \n" << maxRow << endl;
    }

    Array2Xd m = in;
    //        cout<<"m: \n"<< m <<endl; 
    m.row(0) -= maxRow;
    m.row(1) -= maxRow;
    m = m.exp();

    if (nSamples < 20)
    {
        cout << "m: \n" << m << endl;
    }
    out = (w * m.matrix()).array().log();
    out.row(0) += maxRow;
    out.row(1) += maxRow;

    t2 = clock();

    if (nSamples < 20)
    {
        cout << "out \n" << out << endl;
    }

    float diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "TestEigen with " << nSamples << " samples in TIME: " << diff << " sec\n";
}


//void TestEM1()
//{
//    std::cout << "\n###### TestEM\n";
//    SPGM spgm;
//    //    spgm.m_verbose = true;
//
//    int dim = 2;
//    CPTfactor f_a(1, dim);
//    CPTfactor f_a2(1, dim);
//
//    //    CPTfactor f_ba(dim, dim);
//    //    CPTfactor f_ba(dim, dim);
//
//    Vertex a = spgm.AddVNode(0);
//    Vertex a2 = spgm.AddVNode(0);
//
//    //    Vertex prod = spgm.AddProdNode();
//    //    Vertex b = spgm.AddVNode(1);
//    Vertex sum = spgm.AddSumNode();
//    spgm.AddEdge(sum, a);
//    spgm.AddEdge(sum, a2);
//
//    //    Vertex c1 = spgm.AddVNode(2);
//    //    Vertex c2 = spgm.AddVNode(2);
//    //    Vertex d = spgm.AddVNode(3);
//
//    spgm.SetFactor2Root(a, &f_a);
//    spgm.SetFactor2Root(a2, &f_a2);
//
//    //    spgm.SetFactor(b, a, &f_ba);
//    //    spgm.SetFactor(c1, a, &f_ca1);
//    //    spgm.SetFactor(c2, a, &f_ca2);
//    //    spgm.SetFactor(d, c1, &f_dc1);
//    //    spgm.SetFactor(d, c2, &f_dc2);
//
//    //    spgm.AddEdge(a, prod);
//    //    spgm.AddEdge(prod, b);
//    //    spgm.AddEdge(prod, sum);
//    //    spgm.AddEdge(sum, c1);
//    //    spgm.AddEdge(sum, c2);
//    //    spgm.AddEdge(c1, d);
//    //    spgm.AddEdge(c2, d);
//
//    spgm.WriteGraphViz("TestEM1");
//
//    spgm.Finalize();
//
//    SumNode* sn = (SumNode*) spgm.EditNode(sum);
//    sn->SetWeightsRandom();
//
//    Matrix_rowMajor<int> train(3, 1);
//    train.EditRow(0)[0] = 0;
//    train.EditRow(1)[0] = 1;
//    train.EditRow(2)[0] = 1;
//
//
//    //    train.EditRow(1)[0]=1;
//    //    train.EditRow(2)[0]=1;
//
//    //    std::cout<<"FACTORS\n";
//    //    std::cout<<f_a.ToString();
//    //    std::cout<<f_a2.ToString();
//
//    spgm.EM(train, train, 50, 3, 0, true, true);
//
//    cout << spgm.ToString(true);
//
//    for (int i = 0; i < train.nRows(); i++)
//    {
//        cout << "P " << VecToString(train.GetRow(i)) << " = " << VecToStringExponentiated<Real>(spgm.Eval(train.GetRow(i))) << endl;
//    }
//}
//
//void TestEM3()
//{
//    std::cout << "\n###### TestEM\n";
//    SPGM spgm;
//    spgm.m_verbose = false;
//
//    int dim = 2;
//    CPTfactor f_a(1, dim);
//    CPTfactor f_a2(1, dim);
//    CPTfactor f_b(1, dim);
//
//    //    CPTfactor f_ba(dim, dim);
//    //    CPTfactor f_ba(dim, dim);
//
//    Vertex a = spgm.AddVNode(0);
//    Vertex a2 = spgm.AddVNode(0);
//
//    Vertex prod = spgm.AddProdNode();
//    Vertex b = spgm.AddVNode(1);
//    Vertex sum = spgm.AddSumNode();
//    spgm.AddEdge(sum, a);
//    spgm.AddEdge(sum, a2);
//
//    //    Vertex c1 = spgm.AddVNode(2);
//    //    Vertex c2 = spgm.AddVNode(2);
//    //    Vertex d = spgm.AddVNode(3);
//
//    spgm.SetFactor2Root(a, &f_a);
//    spgm.SetFactor2Root(a2, &f_a2);
//    spgm.SetFactor2Root(b, &f_b);
//
//    //    spgm.SetFactor(b, a, &f_ba);
//    //    spgm.SetFactor(c1, a, &f_ca1);
//    //    spgm.SetFactor(c2, a, &f_ca2);
//    //    spgm.SetFactor(d, c1, &f_dc1);
//    //    spgm.SetFactor(d, c2, &f_dc2);
//
//    //    spgm.AddEdge(a, prod);
//    spgm.AddEdge(prod, b);
//    spgm.AddEdge(prod, sum);
//    //    spgm.AddEdge(sum, c1);
//    //    spgm.AddEdge(sum, c2);
//    //    spgm.AddEdge(c1, d);
//    //    spgm.AddEdge(c2, d);
//
//    spgm.WriteGraphViz("TestEM3");
//
//    spgm.Finalize();
//
//    SumNode* sn = (SumNode*) spgm.EditNode(sum);
//    sn->SetWeightsRandom();
//
//    Matrix_rowMajor<int> train(3, 2);
//    train.EditRow(0)[0] = 1;
//    train.EditRow(0)[1] = 0;
//    train.EditRow(1)[0] = 1;
//    train.EditRow(1)[1] = 1;
//    train.EditRow(2)[0] = 1;
//    train.EditRow(2)[1] = 1;
//
//    spgm.EM(train, train, 50, 3, 0, true, true);
//
//    //    std::cout<<"FACTORS\n";
//    //    std::cout<<f_a.ToString();
//    //    std::cout<<f_a2.ToString();
//    //    std::cout<<f_b.ToString();
//
//    cout << spgm.ToString(true);
//
//    for (int i = 0; i < train.nRows(); i++)
//    {
//        cout << "P " << VecToString(train.GetRow(i)) << " = " << VecToStringExponentiated<Real>(spgm.Eval(train.GetRow(i))) << endl;
//    }
//}
//
//void TestEM2()
//{
//    //create a spgm A + (  B   *   (C1 -> D) + (C2 -> D) )
//    SPGM spgm;
//
//
//    int dim = 2;
//    CPTfactor f_a(1, dim);
//    CPTfactor f_ba(dim, dim);
//    CPTfactor f_ca1(dim, dim);
//    CPTfactor f_ca2(dim, dim);
//    CPTfactor f_dc1(dim, dim);
//    CPTfactor f_dc2(dim, dim);
//
//    Vertex a = spgm.AddVNode(0);
//    Vertex prod = spgm.AddProdNode();
//    Vertex b = spgm.AddVNode(1);
//    Vertex sum = spgm.AddSumNode();
//    Vertex c1 = spgm.AddVNode(2);
//    Vertex c2 = spgm.AddVNode(2);
//    Vertex d = spgm.AddVNode(3);
//
//    spgm.SetFactor2Root(a, &f_a);
//    spgm.SetFactor(b, a, &f_ba);
//    spgm.SetFactor(c1, a, &f_ca1);
//    spgm.SetFactor(c2, a, &f_ca2);
//    spgm.SetFactor(d, c1, &f_dc1);
//    spgm.SetFactor(d, c2, &f_dc2);
//
//    spgm.AddEdge(a, prod);
//    spgm.AddEdge(prod, b);
//    spgm.AddEdge(prod, sum);
//    spgm.AddEdge(sum, c1);
//    spgm.AddEdge(sum, c2);
//    spgm.AddEdge(c1, d);
//    spgm.AddEdge(c2, d);
//
//    spgm.WriteGraphViz("TestEM3");
//
//    spgm.Finalize();
//
//    SumNode* sn = (SumNode*) spgm.EditNode(sum);
//    sn->SetWeightsRandom();
//
//    Matrix_rowMajor<int> train(2, 4);
//    train.EditRow(0)[0] = 0;
//    train.EditRow(0)[1] = 1;
//    train.EditRow(0)[2] = 1;
//    train.EditRow(0)[3] = 0;
//
//    train.EditRow(1)[0] = 0;
//    train.EditRow(1)[1] = 1;
//    train.EditRow(1)[2] = 0;
//    train.EditRow(1)[3] = 1;
//
//    //    std::cout<<"FACTORS\n";
//    //    std::cout<<f_a.ToString();
//    //    std::cout<<f_ba.ToString();
//
//    spgm.m_verbose = false;
//    spgm.EM(train, train, 200, 3, 0, true, true);
//
//    cout << spgm.ToString(true);
//
//    for (int i = 0; i < train.nRows(); i++)
//    {
//        cout << "P " << VecToString(train.GetRow(i)) << " = " << VecToStringExponentiated<Real>(spgm.Eval(train.GetRow(i))) << endl;
//    }
//}
//
//void TestEM4()
//{
//    //create a spgm A + (  B   *   (C1 -> D) + (C2 -> D) )
//    SPGM spgm;
//
//
//    int dim = 2;
//    CPTfactor f_a(1, dim);
//    CPTfactor f_ba(dim, dim);
//    CPTfactor f_cb(dim, dim);
//
//    Vertex a = spgm.AddVNode(0);
//    Vertex b = spgm.AddVNode(1);
//    Vertex c = spgm.AddVNode(2);
//
//    spgm.SetFactor2Root(a, &f_a);
//    spgm.SetFactor(b, a, &f_ba);
//    spgm.SetFactor(c, b, &f_cb);
//
//    spgm.AddEdge(a, b);
//    spgm.AddEdge(b, c);
//
//    spgm.WriteGraphViz("TestEM4");
//
//    spgm.Finalize();
//
//    Matrix_rowMajor<int> train(2, 3);
//    train.EditRow(0)[0] = 1;
//    train.EditRow(0)[1] = 1;
//    train.EditRow(0)[2] = 1;
//
//    train.EditRow(1)[0] = 0;
//    train.EditRow(1)[1] = 1;
//    train.EditRow(1)[2] = 0;
//    train.EditRow(1)[3] = 1;
//
//    //    std::cout<<"FACTORS\n";
//    //    std::cout<<f_a.ToString();
//    //    std::cout<<f_ba.ToString();
//
//    spgm.m_verbose = false;
//    spgm.EM(train, train, 200, 3, 0, true, true);
//
//    //    spgm.m_verbose = true;
//    //    
//    //    spgm.Eval(train.GetRow(0));
//    //    spgm.ComputeDerivatives();
//
//    cout << spgm.ToString(true);
//}

void TestChowLiu()
{
    Params p;
    string filename = p.dataFolder + "nltcs.ts.data";
    DataMatrix train;
    train.ReadFromCSV(filename);
    printf("loaded data with %d vars %d samples \n", (int) train.nVars(), (int) train.nData());

    Matrix_rowMajor<Real> MI;
    ProbabilitiesStruct probs;
    std::vector<Real> we = std::vector<Real>();
    ComputeMutualInfo_binary(train, we, MI, probs, 0, true);
    cout << "Mutual Information:\n" << MI.ToString();

    SPGM spgm;
    ChowLiuTree(spgm, MI, probs, true);
    spgm.WriteGraphViz("TestChowLiu");


    Real matlabLL = -6.760055964408658;

    //    spgm.m_verbose = true;
    clock_t t1, t2;
    t1 = clock();
    Real trainLL = spgm.MeanLL(train);
    t2 = clock();
    float diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;

    printf("\n#### RESULTS: \ntrain LL: %.12f, matlab train LL %.12f \nabsolute difference: %.12f\nEVALUATION TIME: %.4f sec", trainLL, matlabLL, abs(trainLL - matlabLL), diff);

    //    spgm.m_verbose = true;
    //    cout << spgm.ToString(true);
    //    spgm.m_verbose = true;
    //    spgm.Eval(train.GetRow(0));

}

void TestCopyConstructor(Params p)
{
    string filename = p.dataFolder + "nltcs.ts.data";
    DataMatrix train;
    train.ReadFromCSV(filename);
    printf("loaded data with %d vars %d samples \n", (int) train.nVars(), (int) train.nData());

    SPGM* spgm = new SPGM();
    LearnMixMaxSpanningTreesSPGM(*spgm, train, p);

    SPGM spgmCopy(*spgm);

    double ll = spgm->MeanLL(train);
    double llCopy = spgmCopy.MeanLL(train);
    printf("\n#### RESULTS: SPGM LL %.12f SPGMcopy LL %.12f difference %.12f\n", ll, llCopy, ll - llCopy);

    delete spgm;
    llCopy = spgmCopy.MeanLL(train);
    printf("\n#### AFTER DELETION: old SPGM LL %.12f SPGMcopy LL %.12f difference %.12f\n", ll, llCopy, ll - llCopy);

    //    printf("\n#### Checking Memory leaks (memory usage should not increase) \n");
    //    for (int i = 0; i < 10; i++)
    //    {
    ////        SPGM other(spgmCopy);
    //        SPGM other = spgmCopy;
    //        spgmCopy = other;
    //        llCopy = other.MeanLL(train); //this does not create leaks
    //        llCopy = spgmCopy.MeanLL(train); //this does
    //    }

}

void TestMixMaxTreesSPGM(Params p)
{
    clock_t t1, t2;
    string filename = p.dataFolder + p.filename + ".ts.data";
    DataMatrix train;
    train.ReadFromCSV(filename);
    printf("loaded data with %d vars %d samples \n", (int) train.nVars(), (int) train.nData());

    SPGM spgm;
    t1 = clock();
    LearnMixMaxSpanningTreesSPGM(spgm, train, p);
    t2 = clock();
    float diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "TestMixMaxTreesSPGM done in " << diff << " sec\n";

    t1 = clock();
    spgm.Eval(train);
    t2 = clock();
    diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "eval done in " << diff << " sec\n";
}

void TestMixMaxTreesSPGM_vsChowLiu(Params p)
{
    string filename = p.dataFolder + "nltcs.ts.data";
    DataMatrix train;
    train.ReadFromCSV(filename);
    printf("loaded data with %d vars %d samples \n", (int) train.nVars(), (int) train.nData());

    string testFile = p.dataFolder + "nltcs.test.data";
    //        string testFile = "testData/test.ts";
    DataMatrix test;
    test.ReadFromCSV(testFile);
    printf("loaded test data with %d vars %d samples \n", (int) test.nVars(), (int) test.nData());

    Matrix_rowMajor<Real> MI;
    ProbabilitiesStruct probs;
    std::vector<Real> we = std::vector<Real>();
    ComputeMutualInfo_binary(train, we, MI, probs, 0, false);

    SPGM clTree;
    ChowLiuTree(clTree, MI, probs, false);
    clTree.WriteGraphViz("TestChowLiu");
    Real trainCL = clTree.MeanLL(train);
    printf("\n#### Chow Liu LL: %.12f\n", trainCL);

    SPGM spgm;
    LearnMixMaxSpanningTreesSPGM(spgm, train, p);
    spgm.WriteGraphViz("TestMixMaxTreesSPGM_vsChowLiu");

    if (p.verbose)
        cout << "############ \n" << spgm.ToString(true);

    clock_t t1, t2;
    spgm.m_verbose = false;

    t1 = clock();
    Real trainSPGM = spgm.MeanLL(train);
    t2 = clock();
    float diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;

    printf("\n#### RESULTS: \nChow Liu LL: %.12f, SPGM LL %.12f \ndifference SPGM - ChowLiu LL : %.12f\nEVALUATION TIME: %.4f sec\n", trainCL, trainSPGM, trainSPGM - trainCL, diff);

    spgm = spgm.EM(train, test, p.nIt_EM, 3, 0, p.em_W, p.em_theta);
    trainSPGM = spgm.MeanLL(test);

    printf("\n#### After EM: \nChow Liu LL: %.12f, SPGM LL %.12f \ndifference SPGM - ChowLiu LL : %.12f\n", trainCL, trainSPGM, trainSPGM - trainCL);


}

//
//void TestRuntime(bool useNltcs, int NmaxPaths, int nTimes)
//{
//    SPGM spgm;
//    std::vector<int> sample;
//
//    if (useNltcs)
//    {
//        string filename = p.dataFolder + "nltcs.ts.data";
//        Matrix_rowMajor<int> train = ReadDataFromCSV(filename);
//        printf("loaded data with %d rows %d cols\n", train.nRows(), train.nCols());
//
//        bool verbose = false;
//        LearnMixMaxSpanningTreesSPGM(spgm, train, NmaxPaths, 0, verbose);
//        spgm.WriteGraphViz("TestRuntime");
//
//        sample = train.GetRow(10);
//    }
//    else
//    {
//        //create a spgm A + (  B   *   (C1 -> D) + (C2 -> D) )
//        spgm.m_verbose = true;
//
//        int dim = 2;
//        CPTfactor* f_a = new CPTfactor(1, dim);
//        CPTfactor* f_ba = new CPTfactor(dim, dim);
//        CPTfactor* f_ca1 = new CPTfactor(dim, dim);
//        CPTfactor* f_ca2 = new CPTfactor(dim, dim);
//        CPTfactor* f_dc1 = new CPTfactor(dim, dim);
//        CPTfactor* f_dc2 = new CPTfactor(dim, dim);
//
//        Vertex a = spgm.AddVNode(0);
//        Vertex prod = spgm.AddProdNode();
//        Vertex b = spgm.AddVNode(1);
//        Vertex sum = spgm.AddSumNode();
//        Vertex c1 = spgm.AddVNode(2);
//        Vertex c2 = spgm.AddVNode(2);
//        Vertex d = spgm.AddVNode(3);
//
//        spgm.SetFactor2Root(a, f_a);
//        spgm.SetFactor(b, a, f_ba);
//        spgm.SetFactor(c1, a, f_ca1);
//        spgm.SetFactor(c2, a, f_ca2);
//        spgm.SetFactor(d, c1, f_dc1);
//        spgm.SetFactor(d, c2, f_dc2);
//
//        spgm.AddEdge(a, prod);
//        spgm.AddEdge(prod, b);
//        spgm.AddEdge(prod, sum);
//        spgm.AddEdge(sum, c1);
//        spgm.AddEdge(sum, c2);
//        spgm.AddEdge(c1, d);
//        spgm.AddEdge(c2, d);
//
//        SumNode* sn = (SumNode*) spgm.EditNode(sum);
//        std::vector<Real> w(2);
//        w[0] = 0.5;
//        w[1] = 0.5;
//        sn->SetWeights(w);
//
//        spgm.Finalize();
//
//        //set sample 
//        sample.push_back(1);
//        sample.push_back(0);
//        sample.push_back(1);
//        sample.push_back(0);
//    }
//    spgm.WriteGraphViz("TestRuntime");
//    clock_t t1, t2;
//    Real trainSPGM;
//
//    spgm.m_verbose = false;
//    t1 = clock();
//    for (int i = 0; i < nTimes; i++)
//        trainSPGM = spgm.Eval(sample);
//    t2 = clock();
//    float diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
//
//    int nMsg = spgm.GetNMessages();
//
//    printf("\n#### TOT MSG EVAL: %d in %.4f sec \n%d messages computed %d times. LL is %f\n\n", nMsg*nTimes, diff, nMsg, nTimes, trainSPGM);
//}

//void TestRuntime2(int nVars, int nSamples)
//{
//    DataMatrix synthData(nVars, nSamples);
//
//    int dim = 2;
//    CPTfactor* f_ba = new CPTfactor(dim, dim);
//    clock_t t1, t2;
//
//    RealMatrix logOutput(2, synthData.nCols());
//
//    t1 = clock();
//    f_ba->Eval(logOutput, synthData, 1, 0);
//    t2 = clock();
//    float diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
//
//    printf("\n#### EVAL with data (%d X %d) in time: %.4f sec \n", nVars, nSamples, diff);
//
//}
//

void TestKnapsack()
{
    vector<Real> vals(3);
    vals[0] = 60;
    vals[1] = 100;
    vals[2] = 120;
    vector<int> costs(3);
    costs[0] = 60;
    costs[1] = 20;
    costs[2] = 30;
    int maxCost = 50;
    vector<bool> chosen;

    float ksv = KnapSack(maxCost, costs, vals, chosen);
    cout << "TestKnapsack()\n";
    cout << "vals: " << VecToString(vals) << endl;
    cout << "costs: " << VecToString(costs) << endl;
    cout << "maxCost: " << maxCost << endl;
    cout << "Knapsack max value: " << ksv << endl;
    cout << "chosen items: " << VecToString(chosen) << endl;
}
//
#include "kmeans-master/kmeans.h"

void TestKMeans(int K)
{
    Params p;
    string filename = p.dataFolder + "nltcs.ts.data";
    DataMatrix train;
    train.ReadFromCSV(filename);
    printf("loaded data with %d vars %d samples \n", (int) train.nVars(), (int) train.nData());

    KMeans kmeans;
    kmeans.run(train, K, 100);
}
//
//void TestSPGM_mixture(int K, int nInsertions, int nEmIters, bool verbose)
//{
//    string trainFile = p.dataFolder + "plants.ts.data";
//    //    string filename = "testData/test.ts";
//    Matrix_rowMajor<int> train = ReadDataFromCSV(trainFile);
//    string testFile = p.dataFolder + "plants.test.data";
//    //    string filename = "testData/test.ts";
//    Matrix_rowMajor<int> valid = ReadDataFromCSV(testFile);
//    printf("loaded data with %d rows %d cols\n", train.nRows(), train.nCols());
//    Real diricheletPrior = 0;
//
//    SPGM_mix spgmMix;
//    LearnMix_MMST_SPGM(spgmMix, train, K, nInsertions, diricheletPrior, verbose);
//
//    Real trainLL = spgmMix.MeanLL(train);
//    printf("\n#### INITAL LL: %.12f\n", trainLL);
//
//    spgmMix.m_verbose = false;
//    spgmMix.EM(train, valid, nEmIters, 3, diricheletPrior, true, true);
//
//    trainLL = spgmMix.MeanLL(train);
//    printf("\n#### FINAL LL %.12f \n", trainLL);
//}
//
//void TestSPGM_mixture2(int K, int nInsertions, int nEmIters, bool verbose)
//{
//    string trainFile = p.dataFolder + "nltcs.ts.data";
//    //    string filename = "testData/test.ts";
//    Matrix_rowMajor<int> train = ReadDataFromCSV(trainFile);
//    string testFile = p.dataFolder + "nltcs.test.data";
//    //    string filename = "testData/test.ts";
//    Matrix_rowMajor<int> valid = ReadDataFromCSV(testFile);
//    printf("loaded data with %d rows %d cols\n", train.nRows(), train.nCols());
//    Real diricheletPrior = 0;
//
//    SPGM_mix spgmMix;
//    LearnMix_MMST_SPGM(spgmMix, train, K, nInsertions, diricheletPrior, verbose);
//
//    Real trainLL = spgmMix.MeanLL(train);
//    printf("\n#### INITAL LL: %.12f\n", trainLL);
//
//    spgmMix.m_verbose = false;
//    spgmMix.EM2(train, valid, nInsertions, nEmIters, 3, diricheletPrior, true, true);
//
//    trainLL = spgmMix.MeanLL(train);
//    printf("\n#### FINAL LL %.12f \n", trainLL);
//}

void TestSPGM_mixture3(Params p)
{
    printf("\n############\nTestSPGM_mixture3(%s, K %d, nInsertions %d, nEmIters1 %d, nEmIters2 %d, verbose %d)\n############\n", p.filename.c_str(), p.K, p.nInsertions, p.EM_iters1, p.EM_iters2, p.verbose);

    string trainFile = p.dataFolder + p.filename + ".ts.data";
    //        string trainFile = "testData/test.ts";
    DataMatrix train;
    train.ReadFromCSV(trainFile);
    printf("loaded %s train with %d vars %d samples \n", p.filename.c_str(), (int) train.nVars(), (int) train.nData());

    string testFile = p.dataFolder + p.filename + ".test.data";
    //        string testFile = "testData/test.ts";
    DataMatrix test;
    test.ReadFromCSV(testFile);
    printf("loaded test data with %d vars %d samples \n", (int) test.nVars(), (int) test.nData());

    SPGM_mix spgmMix(p);
    spgmMix.Init_Kmeans(train, p);

    Real trainLL = spgmMix.MeanLL(test);
    printf("\n#### INITAL test LL: %.12f\n", trainLL);

    spgmMix.m_verbose = false;
    spgmMix.EM_struct(train, test, p);
    spgmMix.EM_params(train, test, p);

    spgmMix.EM_struct(train, test, p);
    spgmMix.EM_params(train, test, p);

    trainLL = spgmMix.MeanLL(test);
    printf("\n#### FINAL test LL %.12f \n", trainLL);
}

void TestSPGM_mixture3_rand(Params p)
{
    string trainFile = p.dataFolder + p.filename + ".ts.data";
    //        string trainFile = "testData/test.ts";
    DataMatrix train;
    train.ReadFromCSV(trainFile);
    printf("loaded %s train with %d vars %d samples \n", p.filename.c_str(), (int) train.nVars(), (int) train.nData());

    string testFile = p.dataFolder + p.filename + ".test.data";
    //        string testFile = "testData/test.ts";
    DataMatrix test;
    test.ReadFromCSV(testFile);
    printf("loaded test data with %d vars %d samples \n", (int) test.nVars(), (int) test.nData());

    SPGM_mix spgmMix(p);
    spgmMix.Init_Random(train, p);

    Real trainLL = spgmMix.MeanLL(test);
    printf("\n#### INITAL test LL: %.12f\n", trainLL);

    spgmMix.m_verbose = false;
    spgmMix.EM_struct(train, test, p);
    spgmMix.EM_params(train, test, p);

    trainLL = spgmMix.MeanLL(test);
    printf("\n#### FINAL test LL %.12f \n", trainLL);
}

void TestSPGM_mixture4(Params p)
{
    string trainFile = p.dataFolder + p.filename + ".ts.data";
    //        string trainFile = "testData/test.ts";
    DataMatrix train;
    train.ReadFromCSV(trainFile);
    printf("loaded train data with %d vars %d samples \n", (int) train.nVars(), (int) train.nData());

    string testFile = p.dataFolder + p.filename + ".test.data";
    //        string testFile = "testData/test.ts";
    DataMatrix test;
    test.ReadFromCSV(testFile);
    printf("loaded test data with %d vars %d samples \n", (int) test.nVars(), (int) test.nData());

    SPGM_mix spgmMix(p);
    spgmMix.Init_Kmeans(train, p);

    Real trainLL = spgmMix.MeanLL(test);
    printf("\n#### INITAL test LL: %.12f\n", trainLL);

    spgmMix.m_verbose = false;
    spgmMix.EM_both(train, test, p);

    spgmMix.EM_params(train, test, p);

    trainLL = spgmMix.MeanLL(test);
    printf("\n#### FINAL test LL %.12f \n", trainLL);
}

void TestSPGM_mixture5(Params p)
{
    string trainFile = p.dataFolder + p.filename + ".ts.data";
    //        string trainFile = "testData/test.ts";
    DataMatrix train;
    train.ReadFromCSV(trainFile);
    printf("loaded train data with %d vars %d samples \n", (int) train.nVars(), (int) train.nData());

    string testFile = p.dataFolder + p.filename + ".test.data";
    //        string testFile = "testData/test.ts";
    DataMatrix test;
    test.ReadFromCSV(testFile);
    printf("loaded test data with %d vars %d samples \n", (int) test.nVars(), (int) test.nData());

    SPGM_mix spgmMix(p);
    spgmMix.Init_Kmeans(train, p);

    Real trainLL = spgmMix.MeanLL(test);
    printf("\n#### INITAL test LL: %.12f\n", trainLL);

    spgmMix.m_verbose = false;

    for (int i = 0; i < p.nIt_EM; i++)
    {
        spgmMix.EM_struct(train, test, p);
        spgmMix.EM_params(train, test, p);
    }

    trainLL = spgmMix.MeanLL(test);
    printf("\n#### FINAL test LL %.12f \n", trainLL);
}