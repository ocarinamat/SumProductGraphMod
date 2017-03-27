/* 
 * File:   DBNBlock.h
 * Author: mdesana
 *
 * Created on June 30, 2014, 4:51 PM
 */

#pragma once

#include "Common.h"
#include "SymmetricSPN.h"
#include "ConvVar.h"
#include "Csi2CsiJunction.h"

/*A DBN block is a structure that constitutes a level in a DBN spn.
 * this is the base class.
 */
class DBNBlock_base
{
protected:

    //this describes each subcase. One case is described by m_OutputElimOrder.
    //The factors are the ones between the input and output vars. 
    //the future order is the order of the output vars (that are eliminated in upper level blocks, not by this one)

    struct SubcaseInfo
    {
        std::vector<ConvVar> m_OutputElimOrder;
        std::vector<int> m_activeIndicesOut; //the indices of the active variables in the output
    };

    std::vector<SubcaseInfo*> m_subcases;

    bool m_structureUpdated;

    const bool m_useTriangulatedDbn; //a parameter

    //    MatrixShards* m_MatrixShards;

    std::vector<ConvVar>* m_inputVars;
    std::vector<ConvVar>* m_outVars;

    DBNBlock_base* m_father;
    FactorCreator* m_factorCreator;
    Csi2CsiJunction* m_csiJunction;

    int m_outW;
    int m_inW;
    int m_level;

    void ComputeEliminOrderInSubnet(int ind);
    void GenerateOutVars();

    bool IsInputLocationActiveInSubnet(int x, SubcaseInfo* subnet)
    {
        assert(!subnet->m_activeIndicesOut.empty());
        int indOut = x;
        std::vector<int>::const_iterator it = std::find(subnet->m_activeIndicesOut.begin(), subnet->m_activeIndicesOut.end(), indOut);
        if (it == subnet->m_activeIndicesOut.end())
            return false;
        else
            return true;
    }

    DBNBlock_base(int inW, int outW, int level);
    DBNBlock_base(const DBNBlock_base& other); //copy constructor
    virtual ~DBNBlock_base();
    

public:
    //adds a subnet, defined by which output variables are active in that case (by position in the output layer) 
    void AddSubnet(const std::vector<Pos2D>& activeOutputPositions); //the coordinates ar in the local grid (starting from 0,0)

    int GetNOutCases() const;
    std::vector<SPN_Var> GetOutVars();
    std::vector<SPN_Var> GetInputVars() const {assert(m_inputVars!=NULL); assert(m_inputVars->size()>0); return ConvertVectorConvVar2SpnVar(*m_inputVars);}  
    const std::vector<ConvVar>* GetOutConvVars();

    std::vector<int> GetActiveIndices(int outCase)
    {
        assert(outCase < GetNOutCases());
        return m_subcases[outCase]->m_activeIndicesOut;
    }

    int GetOutWidth() const
    {
        return m_outW;
    }

    int GetInWidth() const
    {
        return m_inW;
    }

    virtual int GetNofInputCSI() const = 0;

    //Interface

    virtual void Build() = 0;

    virtual void ComputeWeights() = 0; //must be fast
    virtual void IndicatorsSumOut() = 0; //must be fast
    virtual void IndicatorsAssign(const SPN_VarSet &vars); //must be fast 
    virtual void Evaluate() = 0; //must be fast
    virtual int ComputeMaxSPNWidth() = 0;
    virtual std::string ToString() const = 0;

    //related to libDai
    virtual void GetAllFactorsForCase(int csiCase, int outCase, std::vector<FactorFunction*>& f, double* multiplicativeConst) = 0;
    virtual void GetAllEliminationOrderForCase(int csiCase, std::vector<SPN_Var>& order) = 0;
    virtual void GetAllFactors(std::vector<FactorFunction*>& factors, std::vector<FactorFunction*>& csiFactors);
};

//TODO merge all the functionalities of this class into the base class, with the exception of the build method

/* this block is an intermediate block, and has a parent and a child block. 
 * it accepts multiple output cases and multiple input cases (defined by the out cases of the input block)
 */
class DBNBlock : public DBNBlock_base
{
    friend class DBNTopLevelBlock;

protected:
    DBNBlock* m_inputBlock;
    MatrixShards* m_MatrixShards;

public:
    DBNBlock(int inW, int outW, int level); //one shard manager per every input
    DBNBlock(const DBNBlock& other); //copy constructor
    ~DBNBlock();
    //this following 2 functions fill the m_inputVars
    void AddInputBlock(DBNBlock* in); //verifies that this is in exclusion with addition of first level inputs

    //adds a subnet, defined by which output variables are active in that case (by position in the output layer) 
    virtual void Build();

    static void Test();

    void ComputeWeights(); //must be fast
    void IndicatorsSumOut(); //must be fast
//    void IndicatorsAssign(const SPN_VarSet &vars); //must be fast 
    void Evaluate(); //must be fast
    int ComputeMaxSPNWidth();
    std::string ToString() const;

    std::vector<FactorFunction*> GetCSIFactors()
    {
        assert(m_structureUpdated);
        return m_csiJunction->GetCSIFactors();
    }

    CSISumMessage* GetOutCSI(int ind);

    virtual int GetNofInputCSI() const
    {
        assert(m_structureUpdated);
        return m_csiJunction->GetNInputCases();
    }

    //related to libDai
    virtual void GetAllFactorsForCase(int inCase, int outCase, std::vector<FactorFunction*>& f, double* multiplicativeConst);
    virtual void GetAllEliminationOrderForCase(int inCase, std::vector<SPN_Var>& order);
};

/*
 * the first level block is used only in the first level of the dbn spn. 
 * it accepts multiple output cases. 
 * Different from the other blocks, it does not use an input block
 *  */
class DBNFirstLevelBlock : public DBNBlock
{
public:
    DBNFirstLevelBlock(int inW, int outW, int level);
    virtual void Build();
};

/*
 * the top level block is used only in the top level of the dbn spn. 
 * it needs a block in input.
 * Different from the other blocks, its "second layer" is only one variable which can have nClass states 
 *  */
class DBNTopLevelBlock : public DBNBlock_base
{
    DBNBlock* m_inputBlock;
    TableFactorFunction* m_pClass;
    SPN_Var m_classVar;
    SymmetricSPN* m_topSpn;

public:
    DBNTopLevelBlock(int inW, int level, int Nclasses);
    DBNTopLevelBlock(const DBNTopLevelBlock& other);
    virtual ~DBNTopLevelBlock();
    virtual void Build();

    void AddInputBlock(DBNBlock* in);

    virtual void ComputeWeights(); //must be fast
    virtual void IndicatorsSumOut(); //must be fast
    void IndicatorsAssign(const SPN_VarSet &vars); //must be fast 
    virtual void Evaluate(); //must be fast
    virtual int ComputeMaxSPNWidth();
    virtual std::string ToString() const;

    //related to libDai
    void GetAllFactorsForCase(int inCase, int outCase, std::vector<FactorFunction*>& f, double* multiplicativeConst);
    void GetAllEliminationOrderForCase(int inCase, std::vector<SPN_Var>& order);
    virtual void GetAllFactors(std::vector<FactorFunction*>& factors, std::vector<FactorFunction*>& csiFactors);

    std::vector<FactorFunction*> GetCSIFactors()
    {
        assert(m_structureUpdated);
        return m_csiJunction->GetCSIFactors();
    }

    virtual int GetNofInputCSI() const
    {
        assert(m_structureUpdated);
        return m_csiJunction->GetNInputCases();
    }
    double GetRootVal() const;
};