/* 
 * File:   Factor.h
 * Author: mdesana
 *
 * Created on October 14, 2016, 6:33 PM
 */

#pragma once

#include "Common.h"
#include "MatrixContiguous.h"

class CPTfactor
{
public:

    CPTfactor(): m_nOut(0), m_nIn(0){}
    
    CPTfactor(int nOut, int nIn)
    : m_logWeights(nOut, nIn), m_logBeta(nOut, nIn), m_verbose(false), m_nOut(nOut), m_nIn(nIn)
    {
        //         = RealMatrix_rowMajor(nOut, nIn);
        //        m_beta = RealMatrix_rowMajor(nOut, nIn);
        SetRandomWeights();
        ResetBeta();
    }

    int NOut() const
    {
        return m_nOut;
    }

    int NIn() const
    {
        return m_nIn;
    }

    void Eval(RealMatrix& logInput, RealMatrix& logOutput, const DataMatrix& data, const int fromVar, const int toVar) const;
    void Eval(RealMatrix& logOutput, const DataMatrix& data, const int fromVar, const int toVar) const;
    void StackBeta(const RealMatrix& logInput, const RealMatrix& logOutputDerivative, const std::vector<byte>& inputDataRow, const std::vector<Real>& sampleLogProb);
    void PassDerivative_andStackBeta(const RealMatrix& logOutputDerivative, RealMatrix& logInputDerivative, const std::vector<byte>& inputDataRow, const RealMatrix& logInput, const std::vector<Real>& sampleLogProb);

    void SetWeight(int row, int col, Real weight)
    {
        m_logWeights.SetVal(row, col, log(weight));
    }
    
    Real GetLogWeight(int row, int col) const
    {
        return m_logWeights.GetVal(row, col);
    }

    void SetRandomWeights()
    {
        for (int j = 0; j < m_logWeights.nRows(); j++)
        {
            Real sum = 0;
            for (int i = 0; i < m_logWeights.nCols(); i++)
            {
                Real rv = SPN_Rand(RAND_INIT_WEIGHT_MIN, RAND_INIT_WEIGHT_MAX);
                sum += rv;
                m_logWeights.SetVal(j, i, rv);
            }

            Real invSum = 1.f / sum;

            for (int i = 0; i < m_logWeights.nCols(); i++)
            {
                m_logWeights.SetVal(j, i , m_logWeights.GetVal(j, i ) * invSum);
            }
        }
        m_logWeights.ApplyLog();
    }

    void ResetBeta()
    {
        m_logBeta.SetVal(ZERO_LOGVAL);
    }

    std::string ToString()
    {
        std::stringstream s;
        s << "CPTfactor " << this << " (" << m_logWeights.nRows() << " out, " << m_logWeights.nCols() << " in). Weights:\n";
        s << m_logWeights.ToStringExponentiated();
        return s.str();
    }

    void EMupdate(Real uniformPrior);

    static void Test()
    {
        //        int nin = 3;
        //        int nout = 2;
        //        CPTfactor cpt(nout, nin);
        //        cpt.ResetBeta();
        //        cpt.SetRandomWeights();
        //
        //        std::cout << cpt.ToString();
        //
        //        std::vector<Real> logInput(nin);
        //        std::vector<Real> logOutput(nout);
        //
        //        //set values
        //        for (int j = 0; j < cpt.m_logWeights.nRows(); j++)
        //        {
        //            std::vector<Real>& logWrow = cpt.m_logWeights.EditRow(j);
        //            for (int i = 0; i < cpt.m_logWeights.nCols(); i++)
        //                logWrow[i] = log(i + j);
        //        }
        //        for (int i = 0; i < cpt.m_logWeights.nCols(); i++)
        //            logInput[i] = log(i + 1);
        //
        //        std::cout << cpt.ToString();
        //
        //        std::cout << "input\n";
        //        std::cout << VecToStringExponentiated(logInput);
        //
        //        //evaluation
        //        cpt.Eval(logInput, logOutput, -1, 1);
        //        std::cout << "cpt.Eval(logInput,logOutput,1)\n";
        //        std::cout << VecToStringExponentiated(logOutput);
        //
        //        cpt.Eval(logInput, logOutput, -1, -1);
        //        std::cout << "cpt.Eval(logInput,logOutput,-1)\n";
        //        std::cout << VecToStringExponentiated(logOutput);
        //
        //        // derivative
        //        std::vector<Real> logInputDerivative(nin);
        //        std::vector<Real> logOutputDerivative(nout);
        //        for (int i = 0; i < nout; i++)
        //            logOutputDerivative[i] = 0;
        //        for (int i = 0; i < nin; i++)
        //            logInputDerivative[i] = ZERO_LOGVAL;
        //
        //        std::cout << "\nCheck CPTfactor derivative\n";
        //        std::cout << "OutputDerivative\n";
        //        std::cout << VecToStringExponentiated(logOutputDerivative)<<std::endl;
        //        std::cout << "cpt.PassDerivative(logOutputDerivative, logInputDerivative)\n";
        //        cpt.PassDerivative(logOutputDerivative, logInputDerivative, -1);
        //        std::cout << VecToStringExponentiated(logInputDerivative);
        //        
        //        for (int i = 0; i < nout; i++)
        //            logOutputDerivative[i] = ZERO_LOGVAL;
        //        logOutputDerivative[1]=0;
        //        for (int i = 0; i < nin; i++)
        //            logInputDerivative[i] = ZERO_LOGVAL;
        //        cpt.Eval(logInput, logOutput, 2, 1);
        //        std::cout<<"\nCheck derivative with state 1,2\n";
        //        std::cout << "OutputDerivative\n";
        //        std::cout << VecToStringExponentiated(logOutputDerivative)<<std::endl;
        //        std::cout << "Input Derivative\n";
        //        std::cout << VecToStringExponentiated(logInputDerivative)<<std::endl;
        //        std::cout << "cpt.PassDerivative(logOutputDerivative, logInputDerivative)\n";
        //        cpt.PassDerivative(logOutputDerivative, logInputDerivative, 2);
        //        std::cout << VecToStringExponentiated(logInputDerivative);
    }

protected:
    MatrixContiguous<Real> m_logBeta;
    MatrixContiguous<Real> m_logWeights;
    bool m_verbose;
    int m_nOut;
    int m_nIn;

};


