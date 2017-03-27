
#include "Factor.h"

using namespace std;

//todo stack derivative only for right in and out state

void CPTfactor::EMupdate(Real uniformPrior)
{
    bool verbose = m_verbose;

    if (verbose)
        std::cout << "#####EMupdate for following factor \n" << this->ToString() << "with following beta \n" << m_logBeta.ToString() << "\n";

    for (int row = 0; row < m_logBeta.nRows(); row++)
    {
        std::vector<Real> newW = std::vector<Real>(m_logWeights.nCols());

        Real sumTempBeta = 0;

        for (int j = 0; j < m_logWeights.nCols(); j++)
        {
            newW[j] = exp(m_logBeta.GetVal(row, j) + m_logWeights.GetVal(row, j));
            sumTempBeta += newW[j];
        }

        //        if (verbose)
        //        {
        //            std::cout << "factor row " << row << "(beta = " << VecToString(logBetaRow) << ", weights = " << VecToString(logWRow) << ")\n";
        //            for (int j = 0; j < logWRow.size(); j++)
        //            {
        //                printf("        non normalized newW[%d] %f = beta[%d] %f * oldW[%d] %f\n", j, newW[j], j, exp(logBetaRow[j]), j, logWRow[j]);
        //                std::cout << "\n";
        //            }
        //        }

        if (isnan(sumTempBeta))
        {
            std::cerr << "Beta has some inf or nan\n";
        }

        if (sumTempBeta == 0)
        {
            if (verbose)
                std::cout << "no responsibility for this row. continue.\n";
            continue;
        }

        int nw = newW.size();
        Real z = sumTempBeta + uniformPrior;
        Real uOverNw = uniformPrior/nw;
        for (int k = 0; k < newW.size(); k++)
        {
            //            logWRow[k] = log(newW[k] / sumTempBeta);
            m_logWeights.SetVal(row, k, log((newW[k] + uOverNw)/z));
        }

        //DEBUG
        if (verbose)
        {
            std::cout << "Sum beta: " << sumTempBeta;
            std::cout << "New Weights: ";
            for (int k = 0; k < m_logWeights.nCols(); k++)
            {
                printf("%.5f, ", m_logWeights.GetVal(row, k));
            }
            std::cout << "\n";
        }
    }

    ResetBeta();
}

void CPTfactor::StackBeta(const RealMatrix& logInput, const RealMatrix& logOutputDerivative, const std::vector<byte>& inputDataRow, const std::vector<Real>& sampleLogProb)
{
    //for each row, stack derivative
    int N = sampleLogProb.size();
    bool hasInput = !logInput.empty();
    for (int n = 0; n < N; n++)
    {
        byte activeInState = inputDataRow[n];
        Real sampleLL = sampleLogProb[n];
        for (int j = 0; j < m_logBeta.nRows(); j++)
        {
            double logDerivative = logOutputDerivative.GetVal(n, j);
            if (logDerivative != ZERO_LOGVAL)
            {


                double logIn = 0;
                if (hasInput)
                    logIn = logInput.GetVal(n, activeInState);

                //stack the beta for the corresponding weight
                if (logIn != ZERO_LOGVAL)
                {
                    double logBetaQJN = logDerivative + logIn - sampleLL;

                    if (logBetaQJN != ZERO_LOGVAL)
                    {
                        Real& betaRow_activeInState = m_logBeta.EditVal(j, activeInState);
                        if (betaRow_activeInState == ZERO_LOGVAL)
                            betaRow_activeInState = logBetaQJN;
                        else
                            betaRow_activeInState = AddLog(betaRow_activeInState, logBetaQJN);
                        //                        if (m_verbose)
                        //                        {
                        //                            std::printf("    betaRow[%d][%d]+= Derivative[%d] %f*in[%d] %f / sampleP %f)\n", j, activeInState, j, exp(logDerivative), activeInState, exp(logIn), exp(sampleLogProb));
                        //                            std::cout << "    beta is " << m_logBeta.ToString() << endl;
                        //                        }
                    }
#ifndef NDEBUG
                    if (isnan(logBetaQJN))
                    {
                        std::printf("logBetaQJN is nan");
                        std::exit(1);
                    }
#endif
                }
            }

        }
    }
}

//void CPTfactor::StackBeta(const RealMatrix& logInput, const RealMatrix& logOutputDerivative, const std::vector<byte>& inputDataRow, const std::vector<Real>& sampleLogProb)
//{
//    bool hasInput = !logInput.empty();
//    //TODO efficiency reversing for loops
//    Real temp;
//    int N = sampleLogProb.size();
//    //for each row, stack derivative
//    for (int j = 0; j < m_logBeta.nRows(); j++)
//    {
//        const vector<Real>& logDerivative_row_j = logOutputDerivative.GetRow(j);
//        // only for outStateActive or if it is not active TODO
//
//        for (int n = 0; n < N; n++)
//        {
//            Real logDerN = logDerivative_row_j[n];
//            if (logDerN != ZERO_LOGVAL)
//            {
//                byte activeInState = inputDataRow[n];
//                Real logInN = 0;
//                if (hasInput)
//                    logInN = logInput.GetVal(n, activeInState);
//
//                //stack the beta for the corresponding weight
//                if (logInN != ZERO_LOGVAL && logDerN != ZERO_LOGVAL)
//                {
//                    Real logBetaQJN = logDerN + logInN - sampleLogProb[n];
//                    if (logBetaQJN != ZERO_LOGVAL)
//                    {
//                        m_logBeta.AddLogVal(j, activeInState, logBetaQJN);
//                        //                        if (m_verbose)
//                        //                        {
//                        //                            std::printf("    betaRow[%d][%d]+= Derivative[%d] %f*in[%d] %f / sampleP %f)\n", j, activeInState, j, exp(logDerN), activeInState, exp(logInN), exp(sampleLogProb));
//                        //                            std::cout << "    beta is " << m_logBeta.ToString() << endl;
//                        //                        }
//                    }
//#ifndef NDEBUG
//                    if (isnan(logBetaQJN))
//                    {
//
//                        std::printf("logBetaQJN is nan");
//                        std::exit(1);
//                    }
//#endif
//                }
//            }
//        }
//    }
//}

void CPTfactor::PassDerivative_andStackBeta(const RealMatrix& logOutputDerivative, RealMatrix& logInputDerivative, const std::vector<byte>& inputDataRow, const RealMatrix& logInput, const std::vector<Real>& sampleLogProb)
{
    bool hasInput = !logInput.empty();
    int N = inputDataRow.size();

    //for each row, stack derivative
    for (int n = 0; n < N; n++)
    {
        byte activeInState = inputDataRow[n];
        Real currSampleLogProb = sampleLogProb[n];
        for (int j = 0; j < m_logWeights.nRows(); j++)
        {
            //        const std::vector<Real>& logWRow = m_logWeights.GetRow(j);
            //        const std::vector<Real>& logDerivative_row_j = logOutputDerivative.GetRow(j);

            //        if (activeOutState == -1 || activeOutState == j)
            //        {
            Real logD_N = logOutputDerivative.GetVal(n, j);
            if (logD_N != ZERO_LOGVAL)
            {
                //########################## pass derivative

                //update input derivative
                Real l = logD_N + m_logWeights.GetVal(j, activeInState);
                if (m_verbose)
                    std::printf("    l %f = Derivative %f*(weights[%d][%d] %f)\n", exp(l), exp(logD_N), j, activeInState, exp(m_logWeights.GetVal(activeInState, j)));
                Real& inDerN = logInputDerivative.EditVal(n, activeInState);
                if (inDerN == ZERO_LOGVAL)
                    inDerN = l;
                else
                    inDerN = AddLog(inDerN, l);

                //########################## stack beta
                double logIn = 0;
                if (hasInput)
                    logIn = logInput.GetVal(n, activeInState);

                //stack the beta for the corresponding weight
                if (logIn != ZERO_LOGVAL)
                {
                    double logBetaQJN = logD_N + logIn - currSampleLogProb;

                    if (logBetaQJN != ZERO_LOGVAL)
                    {
                        Real& betaRow_activeInState = m_logBeta.EditVal(j, activeInState);
                        if (betaRow_activeInState == ZERO_LOGVAL)
                            betaRow_activeInState = logBetaQJN;
                        else
                            betaRow_activeInState = AddLog(betaRow_activeInState, logBetaQJN);
                        //                        if (m_verbose)
                        //                        {
                        //                            std::printf("    betaRow[%d][%d]+= Derivative[%d] %f*in[%d] %f / sampleP %f)\n", j, activeInState, j, exp(logDerivative), activeInState, exp(logIn), exp(sampleLogProb));
                        //                            std::cout << "    beta is " << m_logBeta.ToString() << endl;
                        //                        }
                    }
#ifndef NDEBUG
                    if (isnan(logBetaQJN))
                    {
                        std::printf("logBetaQJN is nan");
                        std::exit(1);
                    }
#endif
                }
            }
        }
    }
}

//########## 
//########## version where input all at a time (for locality in memory))
//########## 

void CPTfactor::Eval(RealMatrix& logInput, RealMatrix& logOutput, const DataMatrix& data, const int from, const int to) const
{
    assert(logOutput.nCols() == NOut());
    assert(logInput.nCols() == NIn());
    int N = data.nData();
    const ByteMatrix& dataByCol = data.ByCol();

    int in, out;

    if (to == ROOT_VAR)
    {
        for (int n = 0; n < N; n++)
        {
            in = dataByCol.GetVal(from, n);
            logOutput.SetVal(n, 0, m_logWeights.GetVal(0, in) + logInput.GetVal(n, in));
            //todo for debugging, put the rest to NAN
#ifndef NDEBUG
            if (isnan(logOutput.GetVal(n, 0)))
            {
                std::printf("logOutput is nan");
                std::exit(1);
            }
#endif
        }
    }
    else
    {
        for (int n = 0; n < N; n++)
        {
            in = dataByCol.GetVal(from, n);
            out = dataByCol.GetVal(to, n);
            logOutput.SetVal(n, out, m_logWeights.GetVal(out, in) + logInput.GetVal(n, in));
            //todo for debugging, put the rest to NAN
#ifndef NDEBUG
            if (isnan(logOutput.GetVal(n, out)))
            {
                std::printf("logOutput is nan");
                std::exit(1);
            }
#endif
        }
    }
}

void CPTfactor::Eval(RealMatrix& logOutput, const DataMatrix& data, const int from, const int to) const
{
    assert(logOutput.nCols() == NOut());
    int N = data.nData();
    const ByteMatrix& dataByCol = data.ByCol();

    int in, out;
    if (to == ROOT_VAR)
    {
        for (int n = 0; n < N; n++)
        {
            in = dataByCol.GetVal(from, n);
            logOutput.SetVal(n, 0, m_logWeights.GetVal(0, in));
            //todo for debugging, put the rest to NAN
#ifndef NDEBUG
            if (isnan(logOutput.GetVal(n, 0)))
            {
                std::printf("logOutput is nan");
                std::exit(1);
            }
#endif
        }
    }
    else
    {
        for (int n = 0; n < N; n++)
        {
            in = dataByCol.GetVal(from, n);
            out = dataByCol.GetVal(to, n);
            logOutput.SetVal(n, out, m_logWeights.GetVal(out, in));
            //todo for debugging, put the rest to NAN
#ifndef NDEBUG
            if (isnan(logOutput.GetVal(n, out)))
            {
                std::printf("logOutput is nan");
                std::exit(1);
            }
#endif
        }
    }
}

