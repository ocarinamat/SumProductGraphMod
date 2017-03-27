
#include "SPGM_mixture.h"
#include "LearnSPGM.h"
#include "kmeans-master/kmeans.h"

using namespace std;

vector<Real> SPGM_mix::Eval_singleBatch(const DataMatrix& batch, int startIndex)
{
    int Nbatch = batch.nData();
    int N = m_spgmLogVals.nRows();
    assert(m_spgmLogVals.nCols() == m_spgms.size());
    int nInBatch;
    for (int i = 0; i < m_spgms.size(); i++)
    {
        const std::vector<Real>& ll_i = m_spgms[i].Eval(batch);
        nInBatch = 0;
        for (int n = startIndex; n < startIndex + Nbatch; n++)
        {
            assert(n < N);
            m_spgmLogVals.SetVal(n, i, ll_i[nInBatch]);
            nInBatch++;
#ifndef NDEBUG
            assert(!isnan(m_spgmLogVals.GetVal(n, i)));
#endif
        }
    }

    vector<Real> logVals(Nbatch);

    nInBatch = 0;
    for (int n = startIndex; n < startIndex + Nbatch; n++)
    {
        assert(n < N);
        Real logval_n = ZERO_LOGVAL;
        for (int i = 0; i < m_spgms.size(); i++)
        {
            Real logVal_ni = m_spgmLogVals.GetVal(n, i);
            if (logVal_ni != ZERO_LOGVAL)
                logval_n = AddLog(logval_n, logVal_ni + log(m_weights[i]));
        }
        assert(nInBatch < Nbatch);
        logVals[nInBatch] = logval_n;
        nInBatch++;
#ifndef NDEBUG
        assert(!isnan(logval_n));
#endif
    }

    //stack beta of mixture weights
    nInBatch = 0;
    for (int n = startIndex; n < startIndex + Nbatch; n++)
    {
        for (int i = 0; i < m_spgms.size(); i++)
        {
            //stack the beta for this sample
            Real logBetaQJN = m_spgmLogVals.GetVal(n, i) - logVals[nInBatch];
            if (logBetaQJN != ZERO_LOGVAL)
                m_logBeta[i] = AddLog(m_logBeta[i], logBetaQJN);
        }
        nInBatch++;
    }

    return logVals;
}

Real SPGM_mix::Eval(const DataMatrix& data)
{
    return EvalAndDeriv_inBatches(data, false);
}

Real SPGM_mix::EvalAndDeriv_inBatches(const DataMatrix& data, bool doDeriv)
{
    if (m_maxBatchMemSize <= 0)
    {
        cerr << "max batch memory size not set\n";
        exit(1);
    }
    int N = data.nData();
    if (m_spgmLogVals.nRows() != N)
    {
        m_spgmLogVals = RealMatrix(N, m_spgms.size()); //size of the full data 
    }

    //compute memory per sample and determine batch size 
    int memPerSample = 0;
    for (int i = 0; i < m_spgms.size(); i++)
    {
        memPerSample += m_spgms[i].GetMemPerSample();
    }
    int batchSize = floor(double(m_maxBatchMemSize) / double(memPerSample));
    if (batchSize < 1)
        batchSize = 1;
    int nBatches = ceil(double(N) / double(batchSize));

    // eval and deriv for each batch 
    Real meanLL = 0;
    for (int b = 0; b < nBatches; b++)
    {
        int batchInit = b*batchSize;
        const DataMatrix& dataBatch = data.DataSubset(batchInit, batchInit + batchSize);

        //evaluation
        const vector<Real>& ll = SPGM_mix::Eval_singleBatch(dataBatch, batchInit);
        for (int i = 0; i < ll.size(); i++)
            meanLL += ll[i];

        //derivative
        if (doDeriv)
        {
            for (int i = 0; i < m_spgms.size(); i++)
            {
                Real initialLogDerivative = log(m_weights[i]);
                m_spgms[i].ComputeDerivatives(dataBatch, initialLogDerivative, ll); //compute derivatives and beta
            }
        }
    }

    meanLL /= N;
    return meanLL;
}

void SPGM_mix::EMmixtureWeightUpdate(Real uniformPrior)
{
    bool verbose = false;

    std::vector<Real> newW = std::vector<Real>(m_logBeta.size());

    Real sumTempBeta = 0;

    for (int j = 0; j < m_weights.size(); j++)
    {
        newW[j] = exp(m_logBeta[j]) * m_weights[j];
        sumTempBeta += newW[j];
    }

    if (verbose)
    {
        std::cout << "E step in sum node\n";
        for (int j = 0; j < m_weights.size(); j++)
            printf("        non normalized newW[%d] %f = beta[%d] %f * oldW[%d] %f\n", j, newW[j], j, exp(m_logBeta[j]), j, m_weights[j]);
        std::cout << "\n";
    }

    if (isnan(sumTempBeta))
    {
        std::cerr << "Beta has some inf or nan\n";
    }

    ResetMixWeightsBeta();

    if (sumTempBeta == 0)
    {
        assert(false && "no responsibility for this node");
    }

    int nw = newW.size();
    Real z = sumTempBeta + uniformPrior;
    Real uOverNw = uniformPrior / nw;

    for (int k = 0; k < nw; k++)
    {
        m_weights[k] = (newW[k] + uOverNw) / z;
    }

    //DEBUG
    if (verbose)
    {
        std::cout << "Sum beta: " << sumTempBeta;
        std::cout << "New Weights: ";
        for (int k = 0; k < m_weights.size(); k++)
        {
            printf("%.5f, ", m_weights[k]);
        }
        std::cout << "\n";
    }
}

void SPGM_mix::EM_internalSumWeightsOnly(const DataMatrix& train, const DataMatrix& valid, const Params& p)
{
    if (!p.em_W)
        return;

    int maxIters = p.nEmInternalIters;
    int NtriesBeforeEarlyStopping = 0;
    Real uniformWeightPrior = p.uniformWeightPrior;
    bool weightUpdate = p.em_W;
    bool factorsUpdate = false;

    int N = train.nData();
    int worseThanBestValidLLcount = 0;

    if (m_verbose)
        cout << "\n#####learnWithEM\n";

    //compute initial validation LL and save initial params
    if (m_verbose)
        cout << "\n##computing validation LL\n";
    SaveBestParams();
    Real bestValidLL = MeanLL(valid);

    if (m_verbose)
        cout << " initial validation LL " << bestValidLL << endl;

    int bestIter = 0;
    //        bestParams.SaveParams(this);
    int iter;

    Real superOldTrainLL = 1;

    for (iter = 1; iter < maxIters; iter++)
    {
        if (m_verbose)
            cout << "\n##EM iter " << iter << endl;

        //reset the infos for the E step (beta in sum nodes and factors)
        /*#######*/
        for (int i = 0; i < m_spgms.size(); i++)
        {
            m_spgms[i].EM_resetBetas();
            m_logBeta[i] = ZERO_LOGVAL;
        }
        /*#######*/

        //perform forward and backward pass for all samples in training set, saving the info for the E step
        Real oldTrainLL = EvalAndDeriv_inBatches(train, true);

        //then perform the M step for sum nodes and factors

        //weight update for every spn
        //weight update for mixture weights
        /*#######*/
        for (int i = 0; i < m_spgms.size(); i++)
            m_spgms[i].EM_updateStep(uniformWeightPrior, weightUpdate, factorsUpdate);
        //        if (weightUpdate) //we don't do this here
        //        {
        //            EMmixtureWeightUpdate(uniformWeightPrior);
        //        }
        /*#######*/


        if (m_verbose)
            cout << "\n##computing validation LL\n";

        //evaluate validation LL
        Real validLL = MeanLL(valid);

        if (superOldTrainLL != 1 && (superOldTrainLL - oldTrainLL) > 0.00000000000001)
        {
            printf("Warning: training LL decreased from %.12f to %.12f (diff %.13f). This should not happen.\n", superOldTrainLL, oldTrainLL, oldTrainLL - superOldTrainLL);
        }
        superOldTrainLL = oldTrainLL;

        //save best parameters and check convergence
        if (validLL > bestValidLL)
        {
            worseThanBestValidLLcount = 0;
            bestValidLL = validLL;
            SaveBestParams();
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
        printf("      EM_internalSumNode iter %d: oldTrainLL %f, bestValidLL %f, validLL %f\n", iter, oldTrainLL, bestValidLL, validLL);
        if (m_verbose)
            cout << "\n\n";

    }
    if (iter == maxIters)
    {
        printf("EM terminated at iteration %d by reaching maximum iterations\n", iter - 1);
    }

    //set the best params to the SPN 
    LoadBestParams();
    Real currentValidLL = MeanLL(valid);
    if (bestValidLL != currentValidLL)
    {
        printf("ERROR: mismatch in best parameters LL: bestValidLL %f vs currentValidLL %f\n", bestValidLL, currentValidLL);
    }
    printf("best EM validation LL %f was found at iter %d\n\n", currentValidLL, bestIter);
}

Real SPGM_mix::EM_params(const DataMatrix& train, const DataMatrix& valid, const Params& p, const DataMatrix* test /*only for plotting*/)
{
    int maxIters = p.nIt_params;
    int NtriesBeforeEarlyStopping = p.NtriesBeforeEarlyStopping;
    Real uniformWeightPrior = p.uniformWeightPrior;
    bool weightUpdate = p.em_W;
    bool factorsUpdate = p.em_theta;

    int N = train.nData();
    int worseThanBestValidLLcount = 0;

    if (m_verbose)
        cout << "\n#####learnWithEM\n";

    //compute initial validation LL and save initial params
    if (m_verbose)
        cout << "\n##computing validation LL\n";
    SaveBestParams();
    Real bestValidLL = MeanLL(valid);

    if (m_verbose)
        cout << " initial validation LL " << bestValidLL << endl;

    int bestIter = 0;
    //        bestParams.SaveParams(this);
    int iter;

    Real superOldTrainLL = 1;

    for (iter = 1; iter < maxIters; iter++)
    {
        if (m_verbose)
            cout << "\n##EM iter " << iter << endl;

        //reset the infos for the E step (beta in sum nodes and factors)
        /*#######*/
        for (int i = 0; i < m_spgms.size(); i++)
        {
            m_spgms[i].EM_resetBetas();
            m_logBeta[i] = ZERO_LOGVAL;
        }
        /*#######*/

        //perform forward and backward pass for all samples in training set, saving the info for the E step
        Real oldTrainLL = EvalAndDeriv_inBatches(train, true);

        //then perform the M step for sum nodes and factors

        //weight update for every spn
        //weight update for mixture weights
        /*#######*/
        for (int i = 0; i < m_spgms.size(); i++)
            m_spgms[i].EM_updateStep(uniformWeightPrior, weightUpdate, factorsUpdate);
        if (weightUpdate)
        {
            EMmixtureWeightUpdate(uniformWeightPrior);
        }
        /*#######*/


        if (m_verbose)
            cout << "\n##computing validation LL\n";

        //evaluate validation LL
        Real validLL = MeanLL(valid);

        if (superOldTrainLL != 1 && (superOldTrainLL - oldTrainLL) > 0.00000000000001)
        {
            printf("Warning: training LL decreased from %.12f to %.12f (diff %.13f). This should not happen.\n", superOldTrainLL, oldTrainLL, oldTrainLL - superOldTrainLL);
        }
        superOldTrainLL = oldTrainLL;

        //save best parameters and check convergence
        if (validLL > bestValidLL)
        {
            worseThanBestValidLLcount = 0;
            bestValidLL = validLL;
            SaveBestParams();
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
            else
            {
                printf("Validation LL not decreased from optimum. Doing additional trial %d of %d. \n", worseThanBestValidLLcount, NtriesBeforeEarlyStopping);
            }
        }

        if (m_verbose)
            cout << "\n######ITER RESULTS\n";

        //todo remove plotting test LL (it is used only for plotting it)
        if (test != NULL && p.plotTestLL)
        {
            Real testLL = MeanLL(*test);
            printf("  EM_params iter %d: oldTrainLL %f, bestValidLL %f, validLL %f (testLL %f) \n", iter, oldTrainLL, bestValidLL, validLL, testLL);
        }
        else
        {
            printf("  EM_params iter %d: oldTrainLL %f, bestValidLL %f, validLL %f \n", iter, oldTrainLL, bestValidLL, validLL);
        }
        if (m_verbose)
            cout << "\n\n";

    }
    if (iter == maxIters)
    {
        printf("EM terminated at iteration %d by reaching maximum iterations\n", iter - 1);
    }

    //set the best params to the SPN 
    LoadBestParams();
    Real currentValidLL = MeanLL(valid);
    if (bestValidLL != currentValidLL)
    {
        printf("ERROR: mismatch in best parameters LL: bestValidLL %f vs currentValidLL %f\n", bestValidLL, currentValidLL);
    }
    printf("best EM validation LL %f was found at iter %d\n\n", currentValidLL, bestIter);
    return currentValidLL;
}

Real SPGM_mix::EM_struct(const DataMatrix& train, const DataMatrix& valid, const Params& p, const DataMatrix* test /*only for plotting*/)
{
    int N = train.nData();
    int worseThanBestValidLLcount = 0;

    int nInsertions = p.nInsertions;
    int maxIters = p.nIt_struct;
    int NtriesBeforeEarlyStopping = p.NtriesBeforeEarlyStopping;
    Real diricheletPrior = p.diricheletPrior;
    Real lambda = p.lambda;

    if (m_verbose)
        cout << "\n#####learnWithEM\n";

    //compute initial validation LL and save initial params
    if (m_verbose)
        cout << "\n##computing validation LL\n";

    Real bestValidLL = MeanLL(valid);

    if (m_verbose)
        cout << " initial validation LL " << bestValidLL << endl;

    int bestIter = 0;
    SaveBestParams();
    int iter;

    Real superOldTrainLL = 1;

    int K = m_spgms.size();
    Matrix_rowMajor<Real> gamma(K, N);

    for (iter = 1; iter < maxIters; iter++)
    {
        if (m_verbose)
            cout << "\n##EM iter " << iter << endl;

        //compute gammas and oldTrainLL
        Real oldTrainLL = EvalAndDeriv_inBatches(train, false);

        vector<Real> bigGamma(K);
        for (int k = 0; k < K; k++)
            bigGamma[k] = 0;
        for (int n = 0; n < N; n++)
        {
            Real mixVal = 0;
            for (int k = 0; k < m_spgms.size(); k++)
            {
                Real gammaK = m_weights[k] * exp(m_spgmLogVals.GetVal(n, k));
                assert(!isnan(gammaK));
                gamma.EditRow(k)[n] = gammaK;
                mixVal += gammaK;
            }
            for (int k = 0; k < m_spgms.size(); k++)
            {
                if (mixVal == 0)
                    gamma.EditRow(k)[n] = 0;
                else
                    gamma.EditRow(k)[n] /= mixVal;
                assert(!isnan(gamma.GetRow(k)[n]));
                bigGamma[k] += gamma.EditRow(k)[n];
            }
        }
        //################################################

        //then perform the M step for sum nodes and factors
        //compute new weights 
        for (int k = 0; k < K; k++)
            m_weights[k] = bigGamma[k] / N;

        //compute new mixture components
        if (p.fixedSizePerSPGM)
        {
            for (int i = 0; i < K; i++)
            {
                //            delete m_spgms[i];
                m_spgms[i] = SPGM();
                LearnMixMaxSpanningTreesSPGM(m_spgms[i], train, gamma.GetRow(i), p);
            }
        }
        else
        {
            assert(p.takeNbestMi == -1);
            //allow the allocation of edges to be more in one SPGM than in others
            vector<MixMax_auxStruct> auxStruct(K);

            //find candidate paths in all spgms
            for (int i = 0; i < K; i++)
            {
                //            delete m_spgms[i];
                m_spgms[i] = SPGM();
                auxStruct[i] = FindCandidatePaths(m_spgms[i], train, gamma.GetRow(i), p);
            }

            //merge all paths together 
            list<PathInfo> allPaths;
            for (int i = 0; i < K; i++)
            {
                list<PathInfo>& pathsToAdd = auxStruct[i].pathInfosToAdd;
                for (list<PathInfo>::iterator it = pathsToAdd.begin(); it != pathsToAdd.end(); it++)
                {
                    PathInfo path = *it; //todo use ref
                    path.spgmId = i;
                    allPaths.push_back(path);
                }
            }

            //select best paths
            list<PathInfo> selectedPaths = SelectPaths(allPaths, p, true);

            //split paths according to spgm 
            for (list<PathInfo>::iterator it = selectedPaths.begin(); it != selectedPaths.end(); it++)
            {
                PathInfo path = *it; //todo use ref
                int id = path.spgmId;
                auxStruct[id].selectedPathInfos.push_back(path);
            }

            //add paths to proper spgm
            for (int i = 0; i < K; i++)
            {
                AddPaths(m_spgms[i], auxStruct[i], p);
            }
        }

        //add paths to proper spgm
        if (!p.fixedSizePerSPGM)
        {
            cout << "    SPGM sizes: ";
            for (int i = 0; i < K; i++)
                cout << m_spgms[i].GetNMessages() << ", ";
            cout << endl;
        }

        if (p.doParamEmInStructLearn)
        {
            //train ONLY sum node parameters with EM 
            EM_internalSumWeightsOnly(train, valid, p);
        }

        //evaluate validation LL
        if (m_verbose)
            cout << "\n##computing validation LL\n";
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
            SaveBestParams();
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
            else
            {
                printf("Validation LL increased from optimum. Doing additional trial %d of %d. \n", worseThanBestValidLLcount, NtriesBeforeEarlyStopping);
            }
        }

        if (m_verbose)
            cout << "\n######ITER RESULTS\n";

        //todo remove plotting test LL (it is used only for plotting it)
        if (test != NULL && p.plotTestLL)
        {
            Real testLL = MeanLL(*test);
            printf("  EM_struct iter %d: oldTrainLL %f, bestValidLL %f, validLL %f (testLL %f) \n", iter, oldTrainLL, bestValidLL, validLL, testLL);
        }
        else
        {
            printf("  EM_struct iter %d: oldTrainLL %f, bestValidLL %f, validLL %f\n", iter, oldTrainLL, bestValidLL, validLL);
        }

        if (m_verbose)
            cout << "\n\n";

    }
    if (iter == maxIters)
    {
        printf("EM terminated at iteration %d by reaching maximum iterations\n", iter - 1);
    }

    //set the best params to the SPN 
    LoadBestParams();
    Real currentValidLL = MeanLL(valid); //todo avoid last computation
    if (bestValidLL != currentValidLL)
    {
        printf("ERROR: mismatch in best parameters LL: bestValidLL %f vs currentValidLL %f\n", bestValidLL, currentValidLL);
    }
    printf("best EM validation LL %f was found at iter %d\n\n", currentValidLL, bestIter);
    return currentValidLL;
}

void SPGM_mix::EM_both(const DataMatrix& train, const DataMatrix& valid, const Params& p)
{
    int maxIters = p.nIt_params;
    int NtriesBeforeEarlyStopping = p.NtriesBeforeEarlyStopping;
    Real uniformWeightPrior = p.uniformWeightPrior;
    bool weightUpdate = p.em_W;
    bool factorsUpdate = p.em_theta;
    int N = train.nData();
    int worseThanBestValidLLcount = 0;

    if (m_verbose)
        cout << "\n#####learnWithEM\n";

    //compute initial validation LL and save initial params
    if (m_verbose)
        cout << "\n##computing validation LL\n";

    Real bestValidLL = MeanLL(valid);

    if (m_verbose)
        cout << " initial validation LL " << bestValidLL << endl;

    int bestIter = 0;
    SaveBestParams();
    int iter;

    Real superOldTrainLL = 1;

    int K = m_spgms.size();
    Matrix_rowMajor<Real> gamma(K, N);

    for (iter = 1; iter < maxIters; iter++)
    {
        if (m_verbose)
            cout << "\n##EM iter " << iter << endl;

        //compute gammas and oldTrainLL
        Real oldTrainLL = EvalAndDeriv_inBatches(train, false);

        vector<Real> bigGamma(K);
        for (int k = 0; k < K; k++)
            bigGamma[k] = 0;

        for (int n = 0; n < N; n++)
        {
            Real mixVal = 0;
            for (int k = 0; k < m_spgms.size(); k++)
            {
                Real gammaK = m_weights[k] * exp(m_spgmLogVals.GetVal(n, k));
                gamma.EditRow(k)[n] = gammaK;
                mixVal += gammaK;
            }
            for (int k = 0; k < m_spgms.size(); k++)
            {
                gamma.EditRow(k)[n] /= mixVal;
                bigGamma[k] += gamma.EditRow(k)[n];
            }
        }
        //################################################

        //then perform the M step for sum nodes and factors
        //compute new weights 
        for (int k = 0; k < K; k++)
            m_weights[k] = bigGamma[k] / N;

        //compute new mixture components
        for (int i = 0; i < K; i++)
        {
            //            delete m_spgms[i];
            m_spgms[i] = SPGM();
            LearnMixMaxSpanningTreesSPGM(m_spgms[i], train, gamma.GetRow(i), p);
        }

        //##############################################    now perform EM for parameters
        Real validLL = EM_params(train, valid, p);

        //check convergence
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
            SaveBestParams();
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
        printf("  EM super iter %d: oldTrainLL %f, bestValidLL %f, validLL %f\n", iter, oldTrainLL, bestValidLL, validLL);
        if (m_verbose)
            cout << "\n\n";

    }
    if (iter == maxIters)
    {
        printf("EM terminated at iteration %d by reaching maximum iterations\n", iter - 1);
    }

    //set the best params to the SPN 
    LoadBestParams();
    Real currentValidLL = MeanLL(valid);
    if (bestValidLL != currentValidLL)
    {
        printf("ERROR: mismatch in best parameters LL: bestValidLL %f vs currentValidLL %f\n", bestValidLL, currentValidLL);
    }
    printf("best EM validation LL %f was found at iter %d\n\n", currentValidLL, bestIter);
}

void SPGM_mix::Init_Kmeans(const DataMatrix& data, const Params& p)
{
    //run kmeans
    KMeans kmeans;
    Real nMixComponents = p.K;
    const vector<int>& labels = kmeans.run(data, nMixComponents, 50);
    int N = data.nData();
    int dim = data.nVars();
    assert(labels.size() == N);

    assert(m_spgms.empty());

    //    outSpgm = SPGM_mix();
    //    Vertex sumRoot = outSpgm.AddSumNode();
    //    SumNode* sumRootNode = (SumNode*)outSpgm.EditNode(sumRoot);
    //    vector<Real> rootWeights();

    for (int i = 0; i < nMixComponents; i++)
    {
        int NinCluster = 0;
        for (int j = 0; j < labels.size(); j++)
            if (labels[j] == i)
                NinCluster++;

        if (NinCluster == 0)
            continue;

        Real w = (Real) NinCluster / N;
        ByteMatrix bm(NinCluster, dim);
        int currInd = 0;
        for (int j = 0; j < labels.size(); j++)
        {
            if (labels[j] == i)
            {
                bm.EditRow(currInd) = data.ByRow().GetRow(currInd);
                currInd++;
            }
        }
        DataMatrix clusterData;
        clusterData.InitFromMatbyRow(bm);

        SPGM clusterSpgm;
        LearnMixMaxSpanningTreesSPGM(clusterSpgm, clusterData, p);
        AddMixtureElement(w, clusterSpgm);
    }
}

void SPGM_mix::Init_Random(const DataMatrix& data, const Params& p)
{
    int N = data.nData();
    int dim = data.nVars();
    Real nMixComponents = p.K;

    assert(m_spgms.empty());

    //    outSpgm = SPGM_mix();
    //    Vertex sumRoot = outSpgm.AddSumNode();
    //    SumNode* sumRootNode = (SumNode*)outSpgm.EditNode(sumRoot);
    //    vector<Real> rootWeights();

    for (int i = 0; i < nMixComponents; i++)
    {
        std::vector<Real> currWeights(N);
        Real sum = 0;
        for (int n = 0; n < N; n++)
        {
            currWeights[n] = SPN_Rand();
            sum += currWeights[n];
        }
        for (int n = 0; n < N; n++)
            currWeights[n] /= sum;

        SPGM clusterSpgm;
        LearnMixMaxSpanningTreesSPGM(clusterSpgm, data, currWeights, p);
        AddMixtureElement(1.0 / (double) nMixComponents, clusterSpgm);
    }
}