/* 
 * File:   Csi2CsiJunction.h
 * Author: mdesana
 *
 * Created on July 8, 2014, 4:20 PM
 */

#pragma once
#include "Message.h"
#include "ConvVar.h"
#include "SymmetricSPN.h"
#include "FactorCreator.h"

//todo do not use csi if the input cases are only 1

/*
 * a class that joins several symmetric spns into a csi, and prepares several csi nodes to depend on several cases of CSI in output. 
 * This class creates nInCSIXnOutCSI csi symm spns, and there is no sharing between them. so, the shared part but be done outside. 
 */
class Csi2CsiJunction 
{
protected:
    //m_nets are the incoming input networks, one for each input CSI. 
    //we are not responsible for creation and deletion of the messages pointed by m_inMsgsPerCsiCase
    std::vector< std::vector<SumMessageCommon*>* > m_inMsgsPerCsiCase;
    std::vector< std::vector<SymmetricSPN*>* > m_nets; //we handle them internally, create and destroy them 
    std::vector< SPN_Var >  m_nextEliminatedPerOutC; 
    std::vector< std::vector<ProdMessage*>* >  m_prodsPerOutC; //m_prods[i] is a vector containing one per each input CSI. There is one i per each output
    std::vector<CSISumMessage*>   m_outCsiMsgs; //one per each output CSI. The outvars of m_outCsiMsgs[i] are m_varsOutputCases[i]
    std::string m_csiVarPrefix;
    bool m_structureUpdated;
    std::vector<TableFactorFunction*> m_csiFactors; //this is created and deleted by this class
    std::vector<std::vector<ConvVar>* > m_elimOrderPerInCase;
    std::vector<std::vector<ConvVar>* > m_variablesInInMsgPerInCase; //this is a quick fix... ugly
    std::vector<std::vector<ConvVar>* > m_variablesPerOutCase;
    FactorCreator* m_factorCreator;
public:
    Csi2CsiJunction(FactorCreator* f) : m_structureUpdated(false), m_factorCreator(f) {};
    ~Csi2CsiJunction();
    void SetCSIVarPrefix(std::string name){ assert(m_csiVarPrefix.empty()); m_structureUpdated=false; m_csiVarPrefix = name;}
    void AddOutputCase(const std::vector<ConvVar>& outVars, const SPN_Var& nextEliminated); //once inserted, the correct ones to use are selected by the subnets
    void AddInputCase(const std::vector<SumMessageCommon*>& inMsgs, const std::vector<ConvVar>& elimOrder); 
    void AddInputCase(const std::vector<ConvVar>& elimOrder); 
    CSISumMessage* GetOutCSI(int ind) {assert(ind>=0 && ind<m_outCsiMsgs.size()); return m_outCsiMsgs[ind];}
    void Build(); //it also builds the subnets O_o
    int ComputeMaxWidth();
    void ComputeWeights(); //must be fast
    void IndicatorsSumOut(); //must be fast
    void IndicatorsAssign(const SPN_VarSet &vars); //must be fast
    void Evaluate(); //must be fast
    int NOutMessages() const {return m_variablesPerOutCase.size();}
    int GetNInputCases() const {return m_elimOrderPerInCase.size();}
    bool IsBuilt() const {return m_structureUpdated;}
    std::string ToString() const;
    const std::vector<FactorFunction*>& GetFactorsForCase(int inCase, int outCase, double* multiplicativeConstant);
    std::vector<FactorFunction*> GetCSIFactors();
    const std::vector<SPN_Var>& GetEliminationOrder(int inCase) const;
};


