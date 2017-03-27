/* 
 * File:   SPGM_mixture.h
 * Author: mdesana
 *
 * Created on November 11, 2016, 1:59 PM
 */

#pragma once

#include "SPGM.h"

class SPGM_mix
{
    /* this class models a mixture model (convex combination) of SPGMs. 
     Although this could be cast as a single SPGM, merging them would be difficult due to having to replace vertex indices. Therefore we opt for this simple, explicit solution. */

    std::vector<SPGM> m_spgms;
    std::vector<SPGM> m_bestSpgms;
    std::vector<Real> m_bestWeights;
    std::vector<Real> m_weights;
    std::vector<Real> m_logBeta;
    RealMatrix m_spgmLogVals;

    void EMmixtureWeightUpdate(Real uniformPrior);

    void ResetMixWeightsBeta()
    {
        m_logBeta = std::vector<Real>(m_weights.size());
        for (int i = 0; i < m_logBeta.size(); i++)
            m_logBeta[i] = ZERO_LOGVAL;
    }

    void SaveBestParams()
    {
        m_bestWeights = m_weights;
        m_bestSpgms = m_spgms;
    }

    void LoadBestParams()
    {
        m_weights = m_bestWeights;
        m_spgms = m_bestSpgms;
        //        for (int i = 0; i < m_spgms.size(); i++){
        ////            delete m_spgms[i];
        //            m_spgms[i] = m_bestSpgms[i];
        //        }
    }

    Real EvalAndDeriv_inBatches(const DataMatrix& data, bool doDeriv);
    std::vector<Real> Eval_singleBatch(const DataMatrix& batch, int startIndex);
    void EM_internalSumWeightsOnly(const DataMatrix& train, const DataMatrix& valid, const Params& p);
    
public:
    bool m_verbose;
    bool m_debug;
    long int m_maxBatchMemSize;

    SPGM_mix(const Params& p) : m_verbose(false), m_debug(false)
    {
        m_maxBatchMemSize = p.maxBatchMemSize;
    }

    ~SPGM_mix()
    {
        //        for (int i = 0; i < m_spgms.size(); i++)
        //            delete m_spgms[i];
    }

    void AddMixtureElement(Real w, const SPGM& spgm)
    {
        //        assert(spgm != NULL);
        assert(w > 0);
        m_spgms.push_back(spgm);
        m_weights.push_back(w);
        m_logBeta.push_back(NAN);
    }

    Real MeanLL(const DataMatrix& data)
    {
        Real sumw = 0;
        for (int i = 0; i < m_weights.size(); i++)
        {
            sumw += m_weights[i];
            assert(m_weights[i] >= 0);
        }
        if (abs(sumw - 1) > 0.0000000000001)
        {
            assert(false && "non normalized weights");
        }

//        const std::vector<Real>& ll = Eval(data);
//        Real meanLL = 0;
//        for (int i = 0; i < ll.size(); i++)
//            meanLL += ll[i];
//        meanLL /= ll.size();
        return Eval(data);
    }

    Real Eval(const DataMatrix& data);

    Real EM_params(const DataMatrix& training, const DataMatrix& validation, const Params& p, const DataMatrix* test = NULL /*only for plotting*/);
    Real EM_struct(const DataMatrix& training, const DataMatrix& validation, const Params& p, const DataMatrix* test = NULL /*only for plotting*/);
    void EM_both(const DataMatrix& train, const DataMatrix& valid, const Params& p);
    void Init_Kmeans(const DataMatrix& data, const Params& p);
    void Init_Random(const DataMatrix& data, const Params& p);
    static void Test();

    int GetNMessages()
    {
        int n = m_weights.size();
        for (int i = 0; i < m_weights.size(); i++)
            n += m_spgms[0].GetNMessages();
        return n;
    }
};
