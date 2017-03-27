/* 
 * File:   SymmetricSPN.h
 * Author: mdesana
 *
 * Created on May 14, 2014, 5:46 PM
 */

#pragma once

#include "Message.h"

/*this class takes a PGM defined as a set of factors and eliminates them, creating a corresponding SPN.
 * 
 * takes as input a product messages, and eliminates the said nodes creating the corresponding SPN.
 * It creates an output product message which has as indicator the last variable to be eliminated
 * 
 * the elimination order must be specified by the user.
 * Operations:
 *      - create SPN from the graph
 *      - Get the top level message (and ONLY that, the others should be private)
 *      - Assign values to variables
 *      - Evaluate the messages
 *      - compute other thigs such as gradients etc
 *      - use an input an existing message
 * 
 * NOTE: all the functions that are used in the evaluation phase must be fast, the others (used in the construction phase) can also be slower
 * 
 * TODO reuse of messages from an external source?? problem: I might want to distinguish between diferent subnets, torresponding to different CSI
 * todo: (?) is it better not to include the input message indicator into the elimination order? well if included, it is a safety check.. 
 */
//todo rename. This is not a full symmetric SPN but only a part of it (input message and root are missing)
class SymmetricSPN
{
protected:
    
    std::vector<ProdMessage*> m_prods; //ordered
    std::vector<SumMessage*> m_sums; //ordered
    std::vector<SumMessageCommon*> m_unusedSums;
    std::vector<FactorFunction*> m_factors;

    bool m_structureUpdated; //todo usage
    std::vector<FactorFunction*> FindAndEliminateFactors(const SPN_Var& v, std::vector<FactorFunction*>& activeFactors);
    std::vector<SumMessageCommon*> FindAndEliminateMsgByVariables(const SPN_Var& v); //looks in the unused messages which ones depend on a certain variable v
    std::vector<SPN_Var> m_eliminationOrder;
    void Clean();

public:
    SymmetricSPN();
    ~SymmetricSPN();
    void AddFactor(FactorFunction& f);
    void AddInput(SumMessageCommon* sm);
    void SetFactors(const std::vector<FactorFunction*>& factors);
    void AddEliminationInOrder(const SPN_Var& v);
    void SetEliminationOrder(const std::vector<SPN_Var>& varsOrder);
    void Build();
    int ComputeMaxWidth();
    void ComputeWeights(); //must be fast
    void IndicatorsSumOut(); //must be fast
    void IndicatorsAssign(const SPN_VarSet &vars); //must be fast
    void Evaluate(); //must be fast
    const SumMessageCommon* GetRoot() const;
    const std::vector<FactorFunction*>& GetFactors() {return m_factors;}
    double GetRootVal() const {assert(m_structureUpdated); return GetRoot()->GetRootVal();}
    
    const std::vector<SPN_Var>& GetEliminationOrder() const {assert(!m_eliminationOrder.empty()); return m_eliminationOrder;}
    std::vector<SumMessageCommon*> GetOutputMessages() {assert(m_structureUpdated); return m_unusedSums;}
    void GetDepth(); //input depth plus my depth
    bool Eliminates(const SPN_Var& v) const; //check if variable v is eliminated by this spn
    bool IsBuilt() const {return m_structureUpdated;}
    std::string ToString() const;
    void CheckGradient_Weights();
//    void CheckGradient_Activation();
    void ComputeGradient();
};