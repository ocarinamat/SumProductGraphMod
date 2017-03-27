
#include "runs.h"
#include "LearnSPGM.h"
#include <time.h>

using namespace std;

void ReproducePaperResults(Params p)
{
    cout << "#########################################\nREPRODUCING RESULTS OF THE PAPER, WITH BEST PARAMETERS SELECTED WITH CROSS-VALIDATION\n######################################\n\n";
    Real prior;
    p.NtriesBeforeEarlyStopping = 20;

    clock_t t1, t2;
    float diff;
    t1 = clock();

    t2 = clock();
    diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "\n################CLOCK!\nELAPSED TIME " << diff / (3600) << " HOURS####################\n\n ";

    Real tempVal = 1;
    //
    p.filename = "nltcs"; //******
    p.K = tempVal;
    p.nInsertions = tempVal;
    prior = tempVal;
    p.diricheletPrior = prior;
    p.uniformWeightPrior = prior;
    RunOne(p);
    t2 = clock();
    diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "\n################CLOCK!\nELAPSED TIME " << diff / (3600) << " HOURS####################\n\n ";
    //
    p.NtriesBeforeEarlyStopping = 3;
    p.filename = "msnbc"; //******
    p.K = 400;
    p.nInsertions = 5;
    prior = 0.001;
    p.diricheletPrior = prior;
    p.uniformWeightPrior = prior;
    RunOne(p);
    t2 = clock();
    diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "\n################CLOCK!\nELAPSED TIME " << diff / (3600) << " HOURS####################\n\n ";
    //
    p.filename = "kdd";
    p.K = 10;
    p.nInsertions = 200;
    prior = 0.001;
    p.diricheletPrior = prior;
    p.uniformWeightPrior = prior;
    RunOne(p);
    t2 = clock();
    diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "\n################CLOCK!\nELAPSED TIME " << diff / (3600) << " HOURS####################\n\n ";
    //
    p.filename = "plants";
    p.K = 40;
    p.nInsertions = 5000;
    prior = 0.01;
    p.diricheletPrior = prior;
    p.uniformWeightPrior = prior;
    RunOne(p);
    //
    p.filename = "baudio";
    p.K = 20;
    p.nInsertions = 400;
    prior = 0.01;
    p.diricheletPrior = prior;
    p.uniformWeightPrior = prior;
    RunOne(p);
    t2 = clock();
    diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "\n################CLOCK!\nELAPSED TIME " << diff / (3600) << " HOURS####################\n\n ";
    //
    p.filename = "jester";
    p.K = 10;
    p.nInsertions = 1000;
    prior = 0.001;
    p.diricheletPrior = prior;
    p.uniformWeightPrior = prior;
    RunOne(p);
    t2 = clock();
    diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "\n################CLOCK!\nELAPSED TIME " << diff / (3600) << " HOURS####################\n\n ";
    //
    p.filename = "bnetflix";
    p.K = 20;
    p.nInsertions = 400;
    prior = 0.01;
    p.diricheletPrior = prior;
    p.uniformWeightPrior = prior;
    RunOne(p);
    t2 = clock();
    diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "\n################CLOCK!\nELAPSED TIME " << diff / (3600) << " HOURS####################\n\n ";
    //
    p.filename = "accidents"; //******
    p.K = 5;
    p.nInsertions = 1000;
    p.diricheletPrior = 0.001;
    p.uniformWeightPrior = 0.1;
    RunOne(p);
    t2 = clock();
    diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "\n################CLOCK!\nELAPSED TIME " << diff / (3600) << " HOURS####################\n\n ";
    //
    p.filename = "tretail";
    p.K = 5;
    p.nInsertions = 100;
    prior = 0.000000001;
    p.diricheletPrior = prior;
    p.uniformWeightPrior = prior;
    RunOne(p);
    t2 = clock();
    diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "\n################CLOCK!\nELAPSED TIME " << diff / (3600) << " HOURS####################\n\n ";
    //
    p.filename = "pumsb_star"; //******
    p.K = tempVal; //8
    p.nInsertions = tempVal; //5000;
    prior = tempVal; //0.01;
    p.diricheletPrior = prior;
    p.uniformWeightPrior = prior;
    RunOne(p);
    t2 = clock();
    diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "\n################CLOCK!\nELAPSED TIME " << diff / (3600) << " HOURS####################\n\n ";
    //
    p.filename = "dna"; //******
    p.K = 5;
    p.nInsertions = 200;
    prior = 0.001;
    p.diricheletPrior = prior;
    p.uniformWeightPrior = prior;
    RunOne(p);
    t2 = clock();
    diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "\n################CLOCK!\nELAPSED TIME " << diff / (3600) << " HOURS####################\n\n ";

    //
    p.filename = "kosarek";
    p.K = 5;
    p.nInsertions = 1000;
    prior = 0.001;
    p.diricheletPrior = prior;
    p.uniformWeightPrior = prior;
    RunOne(p);
    t2 = clock();
    diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "\n################CLOCK!\nELAPSED TIME " << diff / (3600) << " HOURS####################\n\n ";
    //
    p.filename = "msweb";
    p.K = 5;
    p.nInsertions = 400;
    prior = 0.001;
    p.diricheletPrior = prior;
    p.uniformWeightPrior = prior;
    RunOne(p);
    t2 = clock();
    diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "\n################CLOCK!\nELAPSED TIME " << diff / (3600) << " HOURS####################\n\n ";
    //
    p.filename = "book";
    p.K = 5;
    p.nInsertions = 10;
    prior = 0.001;
    p.diricheletPrior = prior;
    p.uniformWeightPrior = prior;
    RunOne(p);
    t2 = clock();
    diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "\n################CLOCK!\nELAPSED TIME " << diff / (3600) << " HOURS####################\n\n ";
    //
    p.filename = "tmovie";
    p.K = 5;
    p.nInsertions = 5000;
    prior = 0.01;
    p.diricheletPrior = prior;
    p.uniformWeightPrior = prior;
    RunOne(p);
    t2 = clock();
    diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "\n################CLOCK!\nELAPSED TIME " << diff / (3600) << " HOURS####################\n\n ";
    //
    p.filename = "cwebkb";
    p.K = 5;
    p.nInsertions = 5;
    prior = 0.001;
    p.diricheletPrior = prior;
    p.uniformWeightPrior = prior;
    RunOne(p);
    t2 = clock();
    diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "\n################CLOCK!\nELAPSED TIME " << diff / (3600) << " HOURS####################\n\n ";
    //
    p.filename = "cr52";
    p.K = 5;
    p.nInsertions = 400;
    prior = 0.01;
    p.diricheletPrior = prior;
    p.uniformWeightPrior = prior;
    RunOne(p);
    t2 = clock();
    diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "\n################CLOCK!\nELAPSED TIME " << diff / (3600) << " HOURS####################\n\n ";
    //
    p.filename = "c20ng";
    p.K = 5;
    p.nInsertions = 5000;
    prior = 0.01;
    p.diricheletPrior = prior;
    p.uniformWeightPrior = prior;
    RunOne(p);
    t2 = clock();
    diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "\n################CLOCK!\nELAPSED TIME " << diff / (3600) << " HOURS####################\n\n ";
    //
    p.filename = "bbc";
    p.K = 5;
    p.nInsertions = 20;
    prior = 0.1;
    p.diricheletPrior = prior;
    p.uniformWeightPrior = prior;
    RunOne(p);
    t2 = clock();
    diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "\n################CLOCK!\nELAPSED TIME " << diff / (3600) << " HOURS####################\n\n ";
    //
    p.filename = "ad";
    p.K = 200;
    p.nInsertions = 10;
    prior = 0.001;
    p.diricheletPrior = prior;
    p.uniformWeightPrior = prior;
    RunOne(p);
    t2 = clock();
    diff = ((float) t2 - (float) t1) / CLOCKS_PER_SEC;
    cout << "\n################CLOCK!\nELAPSED TIME " << diff / (3600) << " HOURS####################\n\n ";

}

void RunAll(const Params & p)
{
    Params params = p;
    stringstream ss;
    cout << "Running experiments on all datsets\n";
    cout << params.ToString();
    const int Nnames = 20;
    string names[] = {"nltcs", "msnbc", "kdd", "plants", "baudio", "jester", "bnetflix", "accidents", "tretail", "pumsb_star", "dna", "kosarek", "msweb", "book", "tmovie", "cwebkb", "cr52", "c20ng", "bbc", "ad"};

    string dateTime = GetDateTimeString();
    string resultsFile = "RunAll_results_" + dateTime + "_" + p.GetUniqueString() + ".res";

    ss << "Running experiments on all datsets\n";
    ss << params.ToString();
    ss << "\n\nNAME\tTRAIN\tVALID\tTEST\n";
    String2File(ss.str(), resultsFile);
    for (int i = 0; i < Nnames; i++)
    {
        params.filename = names[i];
        Real trainLL;
        Real validLL;
        Real testLL;
        LearnMixEM(params, &trainLL, &validLL, &testLL);
        ss << params.filename << "\t" << trainLL << "\t" << validLL << "\t" << testLL << endl;
        String2File(ss.str(), resultsFile);
    }
    ss << endl;
}

void RunOne(const Params & p)
{
    Params params = p;
    cout << "Running experiment on one datsets\n";
    cout << params.ToString();
    stringstream ss;
    string dateTime = GetDateTimeString();
    string resultsFile = params.filename + "_results_" + dateTime + "_" + p.GetUniqueString() + ".res";

    ss << "Experiment on one dataset\n";
    ss << params.ToString();
    ss << "\n\nNAME\tTRAIN\tVALID\tTEST\n";

    Real trainLL;
    Real validLL;
    Real testLL;
    LearnMixEM(params, &trainLL, &validLL, &testLL);
    ss << params.filename << "\t" << trainLL << "\t" << validLL << "\t" << testLL << endl;
    String2File(ss.str(), resultsFile);
}
