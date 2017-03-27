/* 
 * File:   main1.cpp
 * Author: mdesana
 *
 */

#include "SPGM.h"
#include "tests.h"
#include "runs.h"
#include "LearnSPGM.h"
#include <time.h>

using namespace std;

int main(int argc, char** argv)
{
    Params params;
    params.LoadFromArgs(argc, argv);
    
//    TestMixMaxTreesSPGM(params);

    if (params.reproducePaperResults)
    {
        ReproducePaperResults(params);
    }
    else if (params.runAll)
    {
        RunAll(params);
    }
    else
    {
        RunOne(params);
        
//        cout << params.ToString() << endl;
        
//        
//        Params p = params;
//        
//        p.nInsertions=1;
//        p.em_theta=false;
//        p.diricheletPrior=0;
//        p.uniformWeightPrior=0;
//        
//            string trainFile = p.dataFolder + "nltcs.ts.data";
//    DataMatrix train;
//    train.ReadFromCSV(trainFile);
//    printf("loaded %s train with %d vars %d samples \n", p.filename.c_str(), (int) train.nVars(), (int) train.nData());
//
//    SPGM spgm;
//    LearnMixMaxSpanningTreesSPGM(spgm, train, p);
//    spgm.WriteGraphViz("SPGM_test1");
//    
//    String2File(spgm.ToString(true),"2mt_befEm");
//    spgm = spgm.EM(train, train, 100, 3, 0, true, false); //only train W
//    
//    String2File(spgm.ToString(true),"2mt_aftEm");
        
        
        
        
        

        //    CPTfactor::Test();
        //    SPGMnode::Test();

        //    TestKnapsack();

        //    RealMatrix::Test();
        //        RealMatrix::Test();
        //                SPGM::Test();
        //        TestCopyConstructor(params);
        //            TestChowLiu();
        //        TestMixMaxTreesSPGM(params);

        //    params.lambda = 0;
        //    params.nInsertions = 4;
        //    params.verbose = true;
        //    TestMixMaxTreesSPGM(params);
        //    params.verbose = false;

        //    params.lambda = 0;
        //    params.nInsertions = 120;

        //    string trainFile = params.dataFolder + params.filename + ".ts.data";
        //    Matrix_rowMajor<int> a(4,5);
        //    a.SetVal(1);
        //    a.SetVal(1,2,12);
        //    a.SetVal(0,3,12);
        //    Matrix_rowMajor<int> part = a.GetRowSubset(0,2);
        //    cout<<a.ToString()<<endl;
        //    cout<<part.ToString()<<endl;


        //    params.maxBatchMemSize = 1;
        //    params.nInsertions = 10;
        //    params.K = 2;
        //    params.nIt_struct=10;
        //    params.verbose = true;

        
        
        
        
    }
}

