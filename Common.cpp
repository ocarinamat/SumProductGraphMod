
#include "Common.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>

using namespace std;

Real SPN_Rand(Real fMin, Real fMax)
{
    Real f = (Real) rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

Real SPN_RandSeed()
{
    srand(time(NULL));
}

Real Sigmoid(Real x)
{
    return 1.f / (1.f + exp(-x));
}

void WriteToFile(const std::string& str, const std::string& filename)
{
    std::ofstream file;
    file.open(filename.c_str());
    file << str;
    file.close();
}

// 01 Knapsack with int costs. Returns the maximum value that can be put in a knapsack of capacity maxCost and the configuration for the maximum.

Real KnapSack(int maxCost, std::vector<int>& costs, std::vector<Real>& values, std::vector<bool>& chosenItems)
{
    int n = costs.size();
    assert(costs.size() == values.size());

    Real K[n + 1][maxCost + 1];
    bool keep[n + 1][maxCost + 1];
    for (int i = 0; i <= n; i++)
        for (int w = 0; w <= maxCost; w++)
            keep[i][w] = false;

    // Build table K[][] in bottom up manner
    for (int i = 0; i <= n; i++)
    {
        for (int w = 0; w <= maxCost; w++)
        {
            if (i == 0 || w == 0)
            {
                K[i][w] = 0;
            }
            else if (costs[i - 1] <= w)
            {
                Real valIncluded = values[i - 1] + K[i - 1][w - costs[i - 1]];
                Real valNotIncluded = K[i - 1][w];
                if (valIncluded > valNotIncluded)
                {
                    K[i][w] = valIncluded;
                    keep[i][w] = true;
                }
                else
                {
                    K[i][w] = valNotIncluded;
                    keep[i][w] = false;
                }
            }
            else
            {
                K[i][w] = K[i - 1][w];
                keep[i][w] = false;
            }
        }
    }
    Real maxVal = K[n][maxCost];

    //        cout<<"K: \n[";
    //        for (int r = 0; r < n+1; r++)
    //        {
    //            for (int c = 0; c < maxCost+1; c++)
    //            {
    //                cout<<K[r][c]<<", ";
    //            }
    //            cout<<"]\n";
    //        }
    //        cout<<"\n";
    //        cout<<"keep: \n[";
    //        for (int r = 0; r < n+1; r++)
    //        {
    //            for (int c = 0; c < maxCost+1; c++)
    //            {
    //                cout<<keep[r][c]<<", ";
    //            }
    //            cout<<"]\n";
    //        }
    //        cout<<"\n";

    //compute the chosen items
    chosenItems = std::vector<bool>(n);
    for (int i = 0; i < n; i++)
        chosenItems[i] = false;
    int temp = maxCost;
    for (int i = n; i >= 1; i--)
    {
        if (keep[i][temp] == true)
        {
            temp -= costs[i - 1];
            chosenItems[i - 1] = true;
        }
    }


    //verify that the chosen items have the max val
    Real chosenMaxVal = 0;
    int chosenCost = 0;
    for (int i = 0; i < n; i++)
    {
        if (chosenItems[i])
        {
            chosenMaxVal += values[i];
            chosenCost += costs[i];
        }
    }
    assert(chosenMaxVal == maxVal);
    assert(chosenCost <= maxCost);

    return maxVal;
}

Real AddLog(Real l1, Real l2)
{
    if (l1 > l2)
        return l1 + log(1 + exp(l2 - l1));
    else
        return l2 + log(1 + exp(l1 - l2));
}

std::string SetToString(const std::set<int>& mySet)
{
    std::stringstream ss;
    ss << "{";
    for (std::set<int>::const_iterator it = mySet.begin(); it != mySet.end(); it++)
        ss << *it << ", ";
    ss << "}";
    return ss.str();
}

std::string SetToString(const std::set<Vertex>& mySet)
{
    std::stringstream ss;
    ss << "{";
    for (std::set<Vertex>::const_iterator it = mySet.begin(); it != mySet.end(); it++)
        ss << *it << ", ";
    ss << "}";
    return ss.str();
}

DataMatrix DataMatrix::DataSubset(int startIndex, int stopIndex) const
{
    if (startIndex == 0 && stopIndex == nData())
        return (*this);

    assert(stopIndex >= startIndex);
    if (stopIndex > nData())
        stopIndex = nData();

    int N = stopIndex - startIndex;
    assert(N > 0);

    const ByteMatrix subsetRows = m_byRow.GetRowSubset(startIndex, stopIndex);
    DataMatrix d;
    d.InitFromMatbyRow(subsetRows);
    return d;
}

void DataMatrix::ReadFromCSV(std::string filename)
{
    std::ifstream in(filename.c_str());
    if (!in.is_open())
    {
        cerr << "ERROR: file not found\n";
        Params p;
        p.PrintHelp();
        assert(false && "file not found");
    }

    typedef boost::tokenizer< boost::escaped_list_separator<char> > Tokenizer;

    std::vector< std::string > vec;
    std::string line;

    int nCols = -1;
    int nRows = 0;
    std::list<int> allInts;

    while (std::getline(in, line))
    {
        Tokenizer tok(line);

        vec.assign(tok.begin(), tok.end());
        int nFoundElements = 0;

        for (int i = 0; i < vec.size(); i++)
        {
            try
            {
                boost::trim(vec[i]); //remove trailing spaces
                if (vec[i].empty())
                    continue;
                //                std::cout << "#" << vec[i] << "#";
                int theInt = boost::lexical_cast<short>(vec[i]);
                allInts.push_front(theInt);
                nFoundElements++;
            }
            catch (const boost::bad_lexical_cast &)
            {
                std::cerr << "cannot read int from string: \"" << vec[i] << "\"\n";
                std::cerr << "at row " << nRows << " element " << i << "\n";
                assert(false && "cannot read int");
            }
        }
        //        std::cout << "\n" << std::endl;
        if (nCols == -1)
            nCols = nFoundElements;
        if (nCols != nFoundElements) //all the rows must have same length
        {
            std::cerr << "row " << nRows + 1 << " has " << nFoundElements << " elements " << " but other rows have " << nCols << " elements\n";
            assert(false);
            exit(1);
        }
        nRows++;
    }
    //copy into matrix
    m_byRow = ByteMatrix(nRows, nCols);
    if (nRows == 0 && nCols == 0)
    {
        cerr << "data matrix is empty";
    }
    for (int r = 0; r < nRows; r++)
    {
        for (int c = 0; c < nCols; c++)
        {
            m_byRow.SetVal(r, c, allInts.back() == 1);
            allInts.pop_back();
        }
    }

    assert(allInts.empty());
    m_byCol = m_byRow.Transposed();
}

void String2File(const std::string& s, const std::string filename)
{
    ofstream outFile;
    outFile.open(filename.c_str());
    outFile << s;
    outFile.close();
}

std::string GetDateTimeString()
{
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[200];

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer, 200, "%Y-%m-%d_%I:%M:%S", timeinfo);
    std::string str(buffer);
    return str;
}

//typedef Matrix_rowMajor<Real> RealMatrix;

#include <unistd.h>
#include <getopt.h>

//todo use k means or not?

using namespace std;

void Params::PrintHelp()
{
    cout << "Parameters::\n";
    cout << "   --reproduce (if present, reproduces paper results)"
            "   --K \n"
            "   --nIns \n"
            "   --dir \n"
            "   --wp \n"
            "   --lambda \n"
            "   --init_kmeans(no arg) \n"
            "   --dataFolder \n"
            "   --verbose(no arg) \n"
            "   --emItPar \n"
            "   --emItStruct \n"
            "   --runAll(no arg) \n"
            "   --maxMem \n"
            "   --flexSizes(no arg) \n"
            "   --intEm(no arg) \n"
            "   --takeNbest \n"
            "   --emIt1 \n"
            "   --emIt2 \n"
            "   --plotTestLL \n"
            "   --earlyStopN  \n"
            "   --name : name of dataset in {nltcs, msnbc, kdd, plants, baudio, jester, bnetflix, accidents, tretail, pumsb_star, dna, kosarek, msweb, book, tmovie, cwebkb, cr52, c20ng, bbc, ad} \n";
}

static struct option long_options[] = {
    {"name", required_argument, NULL, 'a'},
    {"emIt1", required_argument, NULL, 'b'},
    {"emIt2", required_argument, NULL, 'c'},
    {"K", required_argument, NULL, 'd'},
    {"nIns", required_argument, NULL, 'e'},
    {"dir", required_argument, NULL, 'f'},
    {"lambda", required_argument, NULL, 'g'},
    {"init_kmeans", no_argument, NULL, 'h'},
    {"dataFolder", required_argument, NULL, 'i'},
    {"verbose", no_argument, NULL, 'j'},
    {"emItPar", required_argument, NULL, 'k'},
    {"emItStruct", required_argument, NULL, 'l'},
    {"runAll", no_argument, NULL, 'm'},
    {"maxMem", required_argument, NULL, 'n'},
    {"wp", required_argument, NULL, 'o'},
    {"flexSizes", no_argument, NULL, 'p'},
    {"intEm", no_argument, NULL, 'q'},
    {"takeNbest", required_argument, NULL, 'r'},
    {"help", no_argument, NULL, 's'},
    {"plotTestLL", no_argument, NULL, 't'},
    {"earlyStopN", required_argument, NULL, 'u'},
    {"reproduce", no_argument, NULL, 'v'},
    {NULL, 0, NULL, 0}
};

void Params::LoadFromArgs(int argc, char** argv)
{
    int ch = 0;

    while ((ch = getopt_long(argc, argv, "abcdefghijklm:", long_options, NULL)) != -1)
    {
        // check to see if a single character or long option came through
        switch (ch)
        {

        case 'a':
            filename = optarg;
            break;
        case 'b':
            EM_iters1 = atoi(optarg);
            break;
        case 'c':
            EM_iters2 = atoi(optarg);
            break;
        case 'd':
            K = atoi(optarg);
            break;
        case 'e':
            nInsertions = atoi(optarg);
            break;
        case 'f':
            diricheletPrior = atof(optarg);
            break;
        case 'g':
            lambda = atof(optarg);
            break;
        case 'h':
            init = Params::init_kmeans;
            break;
        case 'i':
            dataFolder = optarg;
            break;
        case 'j':
            verbose = true;
            break;
        case 'k':
            nIt_params = atoi(optarg);
            break;
        case 'l':
            nIt_struct = atoi(optarg);
            break;
        case 'm':
            runAll = true;
            break;
        case 'n':
            maxBatchMemSize = atol(optarg);
            break;
        case 'o':
            uniformWeightPrior = atof(optarg);
            break;
        case 'p':
            fixedSizePerSPGM = false;
            break;
        case 'q':
            doParamEmInStructLearn = true;
            break;
        case 'r':
            takeNbestMi = atoi(optarg);
            break;
        case 's':
            PrintHelp();
            break;
        case 't':
            plotTestLL = true;
            break;
        case 'u':
            NtriesBeforeEarlyStopping = atoi(optarg);
            break;
        case 'v':
            reproducePaperResults = true;
            break;

        case '?':
            cout << "ERROR: unrecognized parameter\n";
            PrintHelp();
            exit(1);
            break;
        }
    }
}

Params::Params()
{
    //set some default values 

    K = 7;
    diricheletPrior = 0.000000001;
    nInsertions = 140;
    uniformWeightPrior = 0.000000001;
    nEmExternalIters = 200;
    nEmInternalIters = 2;
    EM_iters1 = 200;
    EM_iters2 = 200;
    verbose = false;
    nIt_EM = 200;
    init = Params::init_rand;
    nIt_params = 200;
    nIt_struct = 200;
    NtriesBeforeEarlyStopping = 3;
    filename = "nltcs";
    dataFolder = "/home/mdesana/Downloads/libra-tk-1.1.2c/doc/examples/data/";
    em_theta = true;
    em_W = true;
    lambda = std::numeric_limits<Real>::infinity();
    maxBatchMemSize = 1000000000;
    runAll = false;
    fixedSizePerSPGM = true;
    doParamEmInStructLearn = false;
    maxNAddedPaths = 10000;
    takeNbestMi = -1;
    plotTestLL = false;
    reproducePaperResults = false;
}

string Params::GetUniqueString() const
{
    stringstream ss;
    int initR = init == init_rand;
    ss << setw(4) << "K" << K << "n" << nInsertions << "d" << diricheletPrior << "wp" << uniformWeightPrior << "l" << lambda << "r" << initR << "f" << fixedSizePerSPGM << "w" << em_W << "t" << em_theta << "tn" << takeNbestMi << "nt" << NtriesBeforeEarlyStopping;
    return ss.str();
}

string Params::ToString() const
{
    stringstream ss;
    ss << "#########\nPARAMS:\n";
    ss << "    filename: " << filename << endl;
    ss << "    K: " << K << endl;
    ss << "    nInsertions: " << nInsertions << endl;
    ss << "    diricheletPrior: " << diricheletPrior << endl;
    ss << "    unifWeightPrior: " << uniformWeightPrior << endl;
    ss << "    lambda: " << lambda << endl;
    if (init == init_rand)
        ss << "    init: random\n";
    else if (init == init_kmeans)
        ss << "    init: kmeans\n";

    ss << "    maxBatchMemSize: " << maxBatchMemSize << endl;
    ss << "    em_W: " << em_W << endl;
    ss << "    em_theta: " << em_theta << endl;
    ss << "    nIt_params: " << nIt_params << endl;
    ss << "    nIt_struct: " << nIt_struct << endl;
    ss << "    fixedSizePerSPGM: " << fixedSizePerSPGM << endl;
    ss << "    internalEM_sumW: " << doParamEmInStructLearn << endl;
    ss << "    maxNAddedPaths: " << maxNAddedPaths << endl;
    ss << "    takeNbestMi: " << takeNbestMi << endl;
    ss << "    NtriesBeforeEarlyStopping " << NtriesBeforeEarlyStopping << endl;
    ss << "    reproducePaperResults " << reproducePaperResults << endl;

    ss << "#########";

    return ss.str();

    //    int K;
    //    int nInsertions;
    //    int nEmExternalIters;
    //    int nEmInternalIters;
    //    int EM_iters1;
    //    int EM_iters2;
    //    bool verbose;
    //    int nIt_EM;
    //    int nIt_params;
    //    int nIt_struct;
    //    int NtriesBeforeEarlyStopping;
    //    std::string filename;
    //    std::string dataFolder;
    //    bool em_W;
    //    bool em_theta;
    //    Real diricheletPrior;
}




