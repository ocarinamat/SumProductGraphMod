#include "DBNBlock.h"
#include "FactorCreator.h"

using namespace std;

//################################### base block

DBNBlock_base::DBNBlock_base(int inW, int outW, int level)
: m_useTriangulatedDbn(true)
{
    m_level = level;
    m_inW = inW;
    m_outW = outW;
    m_father = NULL;
    m_structureUpdated = false;
    m_outVars = NULL;
    m_inputVars = NULL;
    m_inputVars = new vector<ConvVar>;
    m_factorCreator = NULL;
    m_csiJunction = NULL;
}

void DBNBlock_base::AddSubnet(const std::vector<Pos2D>& activeOutputPositions)
{
    m_structureUpdated = false;

    //check that all the positions are inside the output area
    for (int i = 0; i < activeOutputPositions.size(); i++)
    {
        assert(activeOutputPositions[i].x < m_outW);
        assert(activeOutputPositions[i].y == 0); //just a convention; the value of y is never used
    }

    SubcaseInfo* sub = new SubcaseInfo;
    std::vector<Pos2D> orderedActivePositions = activeOutputPositions;

    //sort the array so that it corresponds to a good elimination order
    std::sort(orderedActivePositions.begin(), orderedActivePositions.end(), Pos2D::Pos2DSortFunction);

    //now compute the active indices
    vector<int>& actIndicesOut = sub->m_activeIndicesOut;
    assert(actIndicesOut.empty());
    for (int i = 0; i < orderedActivePositions.size(); i++)
    {
        int ind = orderedActivePositions[i].x;
        actIndicesOut.push_back(ind);
    }

    m_subcases.push_back(sub);


#ifdef USE_DEBUG_LOG
    stringstream str;
    str << "added these ordered active positions: ";
    for (int i = 0; i < sub->m_activeIndicesOut.size(); i++)
        str << sub->m_activeIndicesOut[i] << ", ";
    Log::Write(str.str(), LOG_CONV_BLOCK);
#endif
}

DBNBlock_base::~DBNBlock_base()
{
    delete m_inputVars;
    delete m_outVars;
    for (int i = 0; i < m_subcases.size(); i++)
        delete m_subcases[i];
    delete m_factorCreator;
    delete m_csiJunction;
}

void DBNBlock_base::ComputeEliminOrderInSubnet(int i)
{
    assert(!m_inputVars->empty() && !m_outVars->empty()); //call this only after you build the variables

    //the output variables are eliminated by the m_csiJunction. Store them for each network. 
    assert(m_subcases[i]->m_OutputElimOrder.empty());

    const vector<int>& activeInd = m_subcases[i]->m_activeIndicesOut;
    for (int p = 0; p < activeInd.size(); p++)
    {
        int ind = activeInd[p];
        assert(ind >= 0 && ind < m_outVars->size());
        m_subcases[i]->m_OutputElimOrder.push_back(m_outVars->at(ind));
    }
}

void DBNBlock_base::GenerateOutVars()
{
    m_structureUpdated = false;

    //there is one for every input level variable. todo all the inputs must be inserted first. 
    assert(m_outVars == NULL);
    m_outVars = new vector<ConvVar > (m_outW);

    //compute the prefix for variables in this block
    stringstream soutPrefix;
    soutPrefix << "L" << m_level << "_";

#ifdef USE_DEBUG_LOG
    Log::Write("creating all the output level variables for the block", LOG_CONV_BLOCK_VARIABLE_CREATION);
#endif

    assert(!m_inputVars->empty());

    for (int i = 0; i < m_outW; i++)
    {
        ConvVar& outv = m_outVars->at(i);
        outv.m_x = i;
        outv.m_y = 0;

        stringstream sname;
        sname << soutPrefix.str() << outv.m_x;
        SPN_VarManager::PushVar(sname.str(), 2);
        outv.m_var = SPN_VarManager::GetVar(sname.str());

#ifdef USE_DEBUG_LOG
        Log::Write(outv.m_var->GetName(), LOG_CONV_BLOCK_VARIABLE_CREATION);
#endif
    }

    //check that all the variables have been assigned
    vector<ConvVar>& outVars = *(m_outVars);
    for (int i = 0; i < outVars.size(); i++)
    {
        assert(outVars[i].m_var->GetNStates() > 0);
    }
}

int DBNBlock_base::GetNOutCases() const
{
    //    assert(m_structureUpdated);
    assert(m_subcases.size() > 0);
    return m_subcases.size();
}

vector<SPN_Var> DBNBlock_base::GetOutVars()
{
    assert(m_outVars != NULL);
    vector<SPN_Var> outVars;
    vector<ConvVar>& outVarsConv = *m_outVars;
    for (int i = 0; i < outVarsConv.size(); i++)
    {
        outVars.push_back(*(outVarsConv[i].m_var));
    }
    return outVars;
}

const vector<ConvVar>* DBNBlock_base::GetOutConvVars()
{
    return m_outVars;
}

//########################################### DBNBlock

void DBNBlock::AddInputBlock(DBNBlock* in)
{
    assert(in != NULL);
    assert(in->m_father == NULL);
    assert(m_inputBlock == NULL);
    in->m_father = this;

    m_structureUpdated = false;

    m_inputBlock = in;

    const vector<ConvVar>* childOutputVars = in->GetOutConvVars();

    for (int i = 0; i < childOutputVars->size(); i++)
    {
        const ConvVar& cv = childOutputVars->at(i);

        //TODO
        //check for duplicates. It is slow but it's no problem because N^2 with number of variables (small) and only done in building the spn
        for (int j = 0; j < m_inputVars->size(); j++)
            assert(!m_inputVars->at(j).m_var->IsEqual(*cv.m_var));

        //        assert(cv.m_type == ConvVar::poolOrInput);
        m_inputVars->push_back(cv);
    }
}

void DBNTopLevelBlock::AddInputBlock(DBNBlock* in)
{
    assert(in != NULL);
    assert(in->m_father == NULL);
    assert(m_inputBlock == NULL);
    in->m_father = this;
    m_structureUpdated = false;
    m_inputBlock = in;
    //
    //    const vector<ConvVar>* childOutputVars = in->GetOutConvVars();
    //
    //    for (int i = 0; i < childOutputVars->size(); i++)
    //    {
    //        const ConvVar& cv = childOutputVars->at(i);
    //
    //        //TODO
    //        //check for duplicates. It is slow but it's no problem because N^2 with number of variables (small) and only done in building the spn
    //        for (int j = 0; j < m_inputVars->size(); j++)
    //            assert(!m_inputVars->at(j).m_var->IsEqual(*cv.m_var));
    //
    //        //        assert(cv.m_type == ConvVar::poolOrInput);
    //        m_inputVars->push_back(cv);
    //    }
}

DBNBlock::~DBNBlock()
{
    delete m_MatrixShards;
}

DBNBlock::DBNBlock(int inW, int outW, int level)
: DBNBlock_base(inW, outW, level)
{
    m_father = NULL;
    m_inputBlock = NULL;
    m_MatrixShards = new MatrixShards(inW, outW);
    m_inputBlock = NULL;

    //#ifdef USE_DEBUG_LOG
    //    stringstream ss;
    //    ss << "\n######  DBNBlock::DBNBlock:" << this;
    //    ss<<"This is the matrix shard: "<<m_MatrixShards->ToString()<<endl;
    //    Log::Write(ss.str(), LOG_CONV_BLOCK);
    //#endif
}

string DBNBlock::ToString() const
{
    stringstream ss;
    ss << "DBNBlock " << this << endl;
    ss << "There are " << m_csiJunction->NOutMessages() << " output csi cases" << endl;
    ss << "OUTPUT CSI MESSAGES:" << endl;
    for (int i = 0; i < m_csiJunction->NOutMessages(); i++)
    {
        if (m_csiJunction->GetOutCSI(i) != NULL)
        {
            ss << m_csiJunction->GetOutCSI(i) << ": ";
            ss << m_csiJunction->GetOutCSI(i)->ToString() << endl;
        }
        else
            ss << "NULL (no parents in this case)\n";
    }
    return ss.str();
}

void DBNBlock::Build()
{
    stringstream ss;
#ifdef USE_DEBUG_LOG
    Log::Indent();
    ss.str("");
    ss << "\n######  DBNBlock::Build:" << this;
    Log::Write(ss.str(), LOG_CONV_BLOCK);
    Log::Write("USING THIS INPUT BLOCK: ", LOG_CONV_BLOCK);
    ss.str("");
    ss << m_inputBlock << ": ";
    ss << m_inputBlock->ToString() << endl;
    Log::Indent();
    Log::Write(ss.str(), LOG_CONV_BLOCK);
    Log::Dedent();
#endif

    //the computations here do NOT depend on m_father
    assert(!m_subcases.empty());
    //the order of calls here is important
    GenerateOutVars();
    for (int i = 0; i < m_subcases.size(); i++)
    {
        ComputeEliminOrderInSubnet(i);
    }

    assert(!m_inputVars->empty());
    assert(!m_subcases.empty());

    //all the input blocks must have been built, and all their output messages must depend on middle variables
    assert(m_inputBlock != NULL);
    assert(m_inputBlock->m_structureUpdated);

    //    //////////////////////////////////////////////
    //    //in this section we check that the output message of the child block ONLY depends on the output variables of the block TODO move after block build and apply to self
    //    {
    //        bool errorFound = false;
    //        SPN_VarSet outSet;
    //
    //        for (int j = 0; j < m_outVars->size(); j++)
    //            outSet.InsertVariable(*(m_outVars->at(j).m_var));
    //
    //        for (int i = 0; i < m_subcases.size(); i++)
    //        {
    //
    //            //check that the output csi is actually used, and if yes, proceed
    //            if (m_inputBlock->GetOutCSI(i) == NULL)
    //                continue;
    //            const vector<SPN_Var>& inVars = m_inputBlock->GetOutCSI(i)->GetOutVars().GetVars();
    //            for (int v = 0; v < inVars.size(); v++)
    //            {
    //                if (outSet.IndexOf(inVars[v]) < 0)
    //                {
    //                    ss.str("");
    //                    ss << "\n\nERROR: variable " << inVars[v].GetName() << " which is one of the output variables of the input block does not depend on the "
    //                            "set of input variables of this block.\n The input block has address " << m_inputBlock << "\n.This happens in its output case " << i;
    //                    errorFound = true;
    //                    cerr << ss.str();
    //#ifdef USE_DEBUG_LOG
    //                    Log::Write(ss.str(), LOG_CONV_BLOCK);
    //#endif
    //                }
    //            }
    //        }
    //        if (errorFound)
    //        {
    //            cerr << "\n The out variables are:\n" << outSet.ToString_VarsOnly() << "\n\n";
    //            cerr << "\n The in variables are:\n" << outSet.ToString_VarsOnly() << "\n\n";
    //            assert(false);
    //        }
    //    }
    //////////////////////////////////////////////

    //now we create the CSI junction
    assert(m_csiJunction == NULL);

    if (m_useTriangulatedDbn)
    {
        DBN_TriangulatedFactorCreator* newFC = new DBN_TriangulatedFactorCreator(*(m_inputBlock->GetOutConvVars()), *(GetOutConvVars()));
        m_factorCreator = newFC;
    }
    else
    {
        DBN_FactorCreator* newFC = new DBN_FactorCreator(*(m_inputBlock->GetOutConvVars()), *(GetOutConvVars()));
        newFC->m_MatrixShards = m_MatrixShards;
        m_factorCreator = newFC;
    }

    m_csiJunction = new Csi2CsiJunction(m_factorCreator);

    //setup the prefix used in the csi variables in the junction
    stringstream csiPrefix;
    csiPrefix << "L" << m_level << "_CSI";
    m_csiJunction->SetCSIVarPrefix(csiPrefix.str());

#ifdef USE_DEBUG_LOG
    ss.str("");
    ss << "block " << this << " is creating the Csi2CsiJunction " << m_csiJunction;
    Log::Write(ss.str(), LOG_CONV_BLOCK);
#endif
    //first we add the input cases to the csi junction. The input cases are our input messages, and the elimination order is the input elimination order
    for (int i = 0; i < m_subcases.size(); i++)
    {
        vector<SumMessageCommon*> temp;
        temp.push_back(m_inputBlock->GetOutCSI(i));
        if (m_useTriangulatedDbn)
            m_csiJunction->AddInputCase(temp, *(m_inputBlock->GetOutConvVars()));
        else
            m_csiJunction->AddInputCase(temp, m_inputBlock->m_subcases[i]->m_OutputElimOrder);
    }

    int nofcsi = GetNOutCases();
    for (int outc = 0; outc < nofcsi; outc++)
    {
        const SPN_Var& nextElimin = *(m_subcases[outc]->m_OutputElimOrder[0].m_var);
        m_csiJunction->AddOutputCase(m_subcases[outc]->m_OutputElimOrder, nextElimin);
    }

    m_csiJunction->Build();

    m_structureUpdated = true;

#ifdef USE_DEBUG_LOG
    Log::Dedent();
#endif
}

CSISumMessage* DBNBlock::GetOutCSI(int ind)
{
    assert(m_structureUpdated);
    assert(m_csiJunction != NULL);
    assert(m_csiJunction->IsBuilt());
    return m_csiJunction->GetOutCSI(ind);
}

void DBNBlock::ComputeWeights() //must be fast
{
    assert(m_structureUpdated);
    assert(m_csiJunction != NULL);
    m_csiJunction->ComputeWeights();
}

void DBNBlock::IndicatorsSumOut() //must be fast
{
    assert(m_structureUpdated);
    assert(m_csiJunction != NULL);
    m_csiJunction->IndicatorsSumOut();
}

void DBNBlock::Evaluate() //must be fast
{
    assert(m_structureUpdated);
    assert(m_csiJunction != NULL);
    m_csiJunction->Evaluate();
}

void DBNBlock_base::IndicatorsAssign(const SPN_VarSet& vars)
{
    assert(m_structureUpdated);
    assert(m_csiJunction != NULL);
    m_csiJunction->IndicatorsAssign(vars);
}

void DBNTopLevelBlock::IndicatorsAssign(const SPN_VarSet &vars)
{
    DBNBlock_base::IndicatorsAssign(vars);
    m_topSpn->IndicatorsAssign(vars);
}

int DBNBlock::ComputeMaxSPNWidth()
{
    assert(m_structureUpdated);
    assert(m_csiJunction != NULL);
    int max = -1;
    int w = m_csiJunction->ComputeMaxWidth();
    if (w > max)
        max = w;
    return max;
}

void DBNBlock_base::GetAllFactors(std::vector<FactorFunction*>& factors, std::vector<FactorFunction*>& csiFactors)
{
    assert(m_structureUpdated);
    assert(m_factorCreator != NULL);

    const vector<FactorFunction*>& f = m_factorCreator->GetAllFactors();
    factors.insert(factors.end(), f.begin(), f.end());

    const vector<FactorFunction*>& fc = m_csiJunction->GetCSIFactors();
    csiFactors.insert(csiFactors.end(), fc.begin(), fc.end());
}

void DBNTopLevelBlock::GetAllFactors(std::vector<FactorFunction*>& factors, std::vector<FactorFunction*>& csiFactors)
{
    DBNBlock_base::GetAllFactors(factors, csiFactors);
    factors.push_back(m_pClass);
}

void DBNBlock::GetAllFactorsForCase(int inCase, int outCase, std::vector<FactorFunction*>& f, double* multiplicativeConst)
{
    assert(m_structureUpdated);
    const std::vector<FactorFunction*>& caseFactors = m_csiJunction->GetFactorsForCase(inCase, outCase, multiplicativeConst);
    f.insert(f.end(), caseFactors.begin(), caseFactors.end());
}

void DBNBlock::GetAllEliminationOrderForCase(int inCase, std::vector<SPN_Var>& order)
{
    assert(m_structureUpdated);
    const std::vector<SPN_Var>& caseOrder = m_csiJunction->GetEliminationOrder(inCase);
    order.insert(order.end(), caseOrder.begin(), caseOrder.end());
}

//###### DBN first level block

void DBNFirstLevelBlock::Build()
{
    stringstream ss;

#ifdef USE_DEBUG_LOG
    Log::Indent();
    ss.str("");
    ss << "\n######  DBNFirstLevelBlock::Build:" << this;
    Log::Write(ss.str(), LOG_CONV_BLOCK);
#endif

    assert(!m_subcases.empty());
    assert(m_inputBlock == NULL);

    //the order of calls here is important
    GenerateOutVars();
    for (int i = 0; i < m_subcases.size(); i++)
    {
        ComputeEliminOrderInSubnet(i);
    }
    /////////////////////////////////

    assert(!m_inputVars->empty());
    assert(!m_subcases.empty());

    //all the input blocks must have been built, and all their output messages must depend on middle variables

    //now we create the CSI junction
    assert(m_csiJunction == NULL);

    if (m_useTriangulatedDbn)
    {
        DBN_TriangulatedFactorCreator * newFc = new DBN_TriangulatedFactorCreator((*m_inputVars), *(GetOutConvVars()), false);
        m_factorCreator = newFc;
    }
    else
    {
        DBN_FactorCreator* newFc = new DBN_FactorCreator((*m_inputVars), *(GetOutConvVars()));
        newFc->m_MatrixShards = m_MatrixShards;
        m_factorCreator = newFc;
    }


    m_csiJunction = new Csi2CsiJunction(m_factorCreator);

    //setup the prefix used in the csi variables in the junction
    stringstream csiPrefix;
    csiPrefix << "L" << m_level << "_CSI";
    m_csiJunction->SetCSIVarPrefix(csiPrefix.str());

#ifdef USE_DEBUG_LOG
    ss.str("");
    ss << "block " << this << " is creating the Csi2CsiJunction " << m_csiJunction;
    Log::Write(ss.str(), LOG_CONV_BLOCK);
#endif
    //first we add the input cases to the csi junction, that is the full input variables
    m_csiJunction->AddInputCase(*m_inputVars);

    int nofcsi = GetNOutCases();
    for (int outc = 0; outc < nofcsi; outc++)
    {
        const SPN_Var& nextElimin = *(m_subcases[outc]->m_OutputElimOrder[0].m_var);
        m_csiJunction->AddOutputCase(m_subcases[outc]->m_OutputElimOrder, nextElimin);
    }

    m_csiJunction->Build();

    m_structureUpdated = true;

#ifdef USE_DEBUG_LOG
    Log::Dedent();
#endif
}

DBNFirstLevelBlock::DBNFirstLevelBlock(int inW, int outW, int level)
: DBNBlock(inW, outW, level)
{
    m_structureUpdated = false;
    assert(m_inputVars != NULL && m_inputVars->empty()); //and the vectors must be empty

    //create variables in the first level

    for (int x = 0; x < m_inW; x++)
    {
        stringstream str;
        str << "in_" << x;

        ConvVar cvar;
        SPN_VarManager::PushVar(str.str(), 2);
        cvar.m_var = SPN_VarManager::GetVar(str.str());
        cvar.m_x = x;
        cvar.m_y = 0;
        cvar.m_f = 0;
        m_inputVars->push_back(cvar);

#ifdef USE_DEBUG_LOG
        str.str("");
        str << "adding var " << cvar.m_var->GetName() << " with coordinates" << cvar.m_x;
        Log::Write(str.str(), LOG_CONV_BLOCK_VARIABLE_CREATION);
#endif
    }

#ifdef USE_DEBUG_LOG
    Log::Write("", LOG_CONV_BLOCK_VARIABLE_CREATION);
#endif
}

DBNTopLevelBlock::DBNTopLevelBlock(int inW, int level, int Nclasses)
: DBNBlock_base(inW, 1, level)
{
    m_topSpn = NULL;
    m_pClass = NULL;

    //create top var
    SPN_VarManager::PushVar("class", Nclasses);
    m_classVar = *(SPN_VarManager::GetVar("class"));
    m_inputBlock = NULL;
}

DBNTopLevelBlock::~DBNTopLevelBlock()
{
    delete m_topSpn;
    delete m_pClass;
    //    ~DBNBlock_base(); //should be automatically called TODO check
}

void DBNTopLevelBlock::Build()
{
    //all the input blocks must have been built, and all their output messages must depend on middle variables
    assert(m_inputBlock != NULL);
    assert(m_inputBlock->m_structureUpdated);

    stringstream ss;
#ifdef USE_DEBUG_LOG
    Log::Indent();
    ss.str("");
    ss << "\n######  DBNTopLevelBlock::Build:" << this;
    Log::Write(ss.str(), LOG_CONV_BLOCK);
    Log::Write("USING THIS INPUT BLOCK: ", LOG_CONV_BLOCK);
    ss.str("");
    ss << m_inputBlock << ": ";
    ss << m_inputBlock->ToString() << endl;
    Log::Indent();
    Log::Write(ss.str(), LOG_CONV_BLOCK);
    Log::Dedent();
#endif


    //###############################
    //create a csi junction that uses a factor generator, and joins multiple input cases (from the input blocks) to a variable in the output having 
    //nClasses. The latter is summed out later, independently, in a SymmetricSPN

    //now we create the CSI junction
    assert(m_csiJunction == NULL);
    assert(m_factorCreator == NULL);

    ConvVar classConv;
    classConv.m_x = 0;
    classConv.m_var = &m_classVar;
    vector<ConvVar> classVarVec;
    classVarVec.push_back(classConv);

    if (m_useTriangulatedDbn)
    {
        DBN_TriangulatedFactorCreator* fc = new DBN_TriangulatedFactorCreator(*(m_inputBlock->GetOutConvVars()), classVarVec);
        m_factorCreator = fc;
    }
    else
    {
        m_factorCreator = new DBN_TopFactorCreator(*(m_inputBlock->GetOutConvVars()), classConv);
    }
    m_csiJunction = new Csi2CsiJunction(m_factorCreator);

    //generate in a particular way the out vars
    assert(m_outVars == NULL);
    m_outVars = new vector<ConvVar >;
    m_outVars->push_back(classConv);

    //setup the prefix used in the csi variables in the junction
    stringstream csiPrefix;
    csiPrefix << "L" << m_level << "_CSI";
    m_csiJunction->SetCSIVarPrefix(csiPrefix.str());

#ifdef USE_DEBUG_LOG
    ss.str("");
    ss << "block " << this << " is creating the Csi2CsiJunction " << m_csiJunction;
    Log::Write(ss.str(), LOG_CONV_BLOCK);
#endif

    //first we add the input cases to the csi junction. The input cases are our input messages, and the elimination order is the input elimination order
    for (int i = 0; i < m_inputBlock->GetNOutCases(); i++)
    {
        vector<SumMessageCommon*> temp;
        temp.push_back(m_inputBlock->GetOutCSI(i));
        if (m_useTriangulatedDbn)
            m_csiJunction->AddInputCase(temp, *(m_inputBlock->GetOutConvVars()));
        else
            m_csiJunction->AddInputCase(temp, m_inputBlock->m_subcases[i]->m_OutputElimOrder);
    }

    m_csiJunction->AddOutputCase(classVarVec, m_classVar);

    m_csiJunction->Build();

    //now create the symmetric SPN summing out the class variable
    assert(m_topSpn == NULL);
    m_topSpn = new SymmetricSPN;
    assert(m_pClass == NULL);

    SPN_VarSet v;
    v.InsertVariable(m_classVar);
    m_pClass = new TableFactorFunction(v);
    m_pClass->SetRandom();
    m_topSpn->AddFactor(*m_pClass);

    assert(m_csiJunction->NOutMessages() == 1);
    m_topSpn->AddInput(m_csiJunction->GetOutCSI(0));

    vector<SPN_Var> topEliminOrder;
    topEliminOrder.push_back(m_classVar);
    m_topSpn->SetEliminationOrder(topEliminOrder);

    //now build the top level SPN
    m_topSpn->Build();

    m_structureUpdated = true;
}

void DBNTopLevelBlock::ComputeWeights()
{
    assert(m_structureUpdated);
    assert(m_csiJunction != NULL);
    m_csiJunction->ComputeWeights();
    m_topSpn->ComputeWeights();
}

void DBNTopLevelBlock::IndicatorsSumOut()
{
    assert(m_structureUpdated);
    assert(m_csiJunction != NULL);
    m_csiJunction->IndicatorsSumOut();
    m_topSpn->IndicatorsSumOut();
}

void DBNTopLevelBlock::Evaluate()
{
    assert(m_structureUpdated);
    assert(m_csiJunction != NULL);
    m_csiJunction->Evaluate();
    m_topSpn->Evaluate();
}

int DBNTopLevelBlock::ComputeMaxSPNWidth()
{
    assert(m_structureUpdated);
    assert(m_csiJunction != NULL);
    int max = -1;
    int w = m_csiJunction->ComputeMaxWidth();
    if (w > max)
        max = w;
    w = m_topSpn->ComputeMaxWidth();
    if (w > max)
        max = w;
    return max;
}

std::string DBNTopLevelBlock::ToString() const
{
    stringstream s;
    s << "top level block with class variable " << m_classVar.GetName(); //TODO
    return s.str();
}

void DBNTopLevelBlock::GetAllFactorsForCase(int inCase, int outCase, std::vector<FactorFunction*>& f, double* multiplicativeConst)
{
    assert(outCase == 0 && "there is only one output case in the top layer");
    assert(m_structureUpdated);
    const std::vector<FactorFunction*>& caseFactors = m_csiJunction->GetFactorsForCase(inCase, outCase, multiplicativeConst);
    f.insert(f.end(), caseFactors.begin(), caseFactors.end());

    //we also add the class unary factor
    f.push_back(m_pClass);
}

void DBNTopLevelBlock::GetAllEliminationOrderForCase(int inCase, std::vector<SPN_Var>& order)
{
    assert(m_structureUpdated);
    const std::vector<SPN_Var>& caseOrder = m_csiJunction->GetEliminationOrder(inCase);
    order.insert(order.end(), caseOrder.begin(), caseOrder.end());
    order.push_back(m_classVar);
}

double DBNTopLevelBlock::GetRootVal() const
{
    assert(m_structureUpdated);
    return m_topSpn->GetRootVal();
}