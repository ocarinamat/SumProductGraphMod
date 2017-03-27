#include "Csi2CsiJunction.h"
#include "FactorCreator.h"

using namespace std;

int Csi2CsiJunction::ComputeMaxWidth()
{
    assert(m_structureUpdated);
    int max = -1;
    for (int i = 0; i < m_nets.size(); i++)
    {
        for (int j = 0; j < m_nets[i]->size(); j++)
        {
            int temp = m_nets[i]->at(j)->ComputeMaxWidth();
            if (temp > max)
                max = temp;
        }
    }
    //    for (int i = 0; i < m_prods.size(); i++)
    //    {
    //        int temp = m_prods[i]->ComputeWidth();
    //        if (temp > max)
    //            max = temp;
    //    }
    return max;
}

const std::vector<FactorFunction*>& Csi2CsiJunction::GetFactorsForCase(int inCase, int outCase, double* multiplicativeConstant)
{
    assert(m_structureUpdated);
    int nInCases = GetNInputCases();
    assert(inCase < nInCases);
    assert(outCase < NOutMessages());

    //the multiplicative constant is the value of the csi factor associated to a particular out case, when the csi variable has state inCase
    assert(multiplicativeConstant != NULL);
    assert(outCase < m_csiFactors.size());
    assert(m_csiFactors[outCase] != NULL);
    vector<int> assignment(1, inCase);
    *multiplicativeConstant = m_csiFactors[outCase]->Compute(assignment);

    return m_nets[outCase]->at(inCase)->GetFactors();
}

const std::vector<SPN_Var>& Csi2CsiJunction::GetEliminationOrder(int inCase) const
{
    assert(m_structureUpdated);
    int nInCases = GetNInputCases();
    assert(inCase < nInCases);
    // I take the output case 0 network because they all have the same elimination order (that only depends on the input case - see Csi2CsiJunction::AddInputCase)
    return m_nets[0]->at(inCase)->GetEliminationOrder();
}

string Csi2CsiJunction::ToString() const
{
    stringstream s;
    s << "\n###########  Csi2CsiJunction INFO: ID " << this << "\n";
    s << "Nets:\n";
    for (int o = 0; o < m_nets.size(); o++)
    {
        for (int i = 0; i < m_nets.size(); i++)
        {
            s << "net for output " << o << " input " << i << " ID " << m_nets[o]->at(i) << endl << m_nets[o]->at(i)->ToString() << endl;
        }
    }
    s << "\n END OF Csi2CsiJunction INFO #######";
    return s.str();
}

std::vector<FactorFunction*> Csi2CsiJunction::GetCSIFactors()
{
    assert(m_structureUpdated);
    std::vector<FactorFunction*> f;
    for (int i = 0; i < m_csiFactors.size(); i++)
        f.push_back(m_csiFactors[i]);
    return f;
}

Csi2CsiJunction::~Csi2CsiJunction()
{
    for (int i = 0; i < m_outCsiMsgs.size(); i++)
        delete m_outCsiMsgs[i];

    for (int i = 0; i < m_nets.size(); i++)
    {
        for (int j = 0; j < m_nets[i]->size(); j++)
        {
            //            if (m_nets[i]->at(j) != NULL)
            //            {
            //                vector<FactorFunction*> f = (m_nets[i]->at(j)->GetFactors());
            //                for (int k = 0; k < f.size(); k++)
            //                    delete f[k];
            //            }
            delete m_nets[i]->at(j);
        }
        delete m_nets[i];
    }

    for (int i = 0; i < m_prodsPerOutC.size(); i++)
    {
        if (m_prodsPerOutC[i] != NULL)
        {
            for (int j = 0; j < m_prodsPerOutC[i]->size(); j++)
                delete m_prodsPerOutC[i]->at(j);
            delete m_prodsPerOutC[i];
        }
    }

    for (int i = 0; i < m_elimOrderPerInCase.size(); i++)
    {
        delete m_elimOrderPerInCase[i];
    }

    for (int i = 0; i < m_csiFactors.size(); i++)
    {
        delete m_csiFactors[i];
    }

    for (int i = 0; i < m_inMsgsPerCsiCase.size(); i++)
        delete m_inMsgsPerCsiCase[i];

    for (int i = 0; i < m_variablesPerOutCase.size(); i++)
        delete m_variablesPerOutCase[i];
}

void Csi2CsiJunction::AddOutputCase(const std::vector<ConvVar>& outVars, const SPN_Var& nextEliminated)
{
    const vector<SPN_Var>& vars = ConvertVectorConvVar2SpnVar(outVars);
    SPN_VarSet vs(vars);

#ifdef USE_DEBUG_LOG
    Log::Indent();
    stringstream s;
    s << "\n####################### Csi2CsiJunction::AddOutputCase \n vars: ";
    s << vs.ToString_VarsOnly();
    Log::Write("NEXT eliminated var:", LOG_CSI_JUNCTION);
    Log::Write(nextEliminated.GetName(), LOG_CSI_JUNCTION);
    s << "\n####################### END AddOutputCase\n";
    Log::Write(s.str(), LOG_CSI_JUNCTION);

    Log::Dedent();
#endif

    assert(!outVars.empty());
    m_structureUpdated = false;
    m_nextEliminatedPerOutC.push_back(nextEliminated);
    vector<ConvVar>* outv = new vector<ConvVar>;
    *outv = outVars;
    m_variablesPerOutCase.push_back(outv);

    assert(vs.IndexOf(nextEliminated) >= 0 && "the next eliminated must be part of the output variables");

    //    //check that the next eliminated is in some of the junction factors //todo? move this in the factor generating function?
    //    bool found = false;
    //    for (int i = 0; i < junctionFactors.size(); i++)
    //    {
    //        int ind = junctionFactors[i]->GetVarSet().IndexOf(nextEliminated);
    //        if (ind >= 0)
    //        {
    //            found = true;
    //            break;
    //        }
    //    }
    //    if (!found)
    //    {
    //        stringstream msg;
    //        msg << "\nERROR: the next eliminated variable " << nextEliminated.GetName() << "is not in the variables of the juction factors\n\n";
    //        cerr << msg.str();
    //        assert(false);
    //    }
}

void Csi2CsiJunction::AddInputCase(const std::vector<SumMessageCommon*>& inMsgs, const std::vector<ConvVar>& elimOrder)
{
#ifdef USE_DEBUG_LOG
    Log::Indent();
    stringstream s;
    s << "\n####################### Csi2CsiJunction::AddInputCase \n Eliminated vars in order: ";
    s << "Eliminated vars in order: ";
    for (int i = 0; i < elimOrder.size(); i++)
        s << elimOrder[i].m_var->GetName() << " ";
    s << "\n####################### END AddInputCase\n";
    Log::Write(s.str(), LOG_CSI_JUNCTION);
    Log::Dedent();
#endif

    m_structureUpdated = false;
    std::vector<SumMessageCommon*>* copyInMsg = new std::vector<SumMessageCommon*>;
    *copyInMsg = inMsgs;
    m_inMsgsPerCsiCase.push_back(copyInMsg); //inMsgs can be null (in which case no input messages are provided
    std::vector<ConvVar>* v = new std::vector<ConvVar>;
    *v = elimOrder;
    m_elimOrderPerInCase.push_back(v);

    //find the variables that are in elimOrder AND in inMsgs
    std::vector<ConvVar>* w = new std::vector<ConvVar>;
    SPN_VarSet inMsgVars;
    for (int i = 0; i < inMsgs.size(); i++)
        inMsgVars.InsertVariables(inMsgs[i]->GetOutVars().GetVars());
    for (int i = 0; i < elimOrder.size(); i++)
        if (inMsgVars.IndexOf(elimOrder[i].GetVar()) >= 0)
            w->push_back(elimOrder[i]);
    m_variablesInInMsgPerInCase.push_back(w);
}

void Csi2CsiJunction::AddInputCase(const std::vector<ConvVar>& elimOrder)
{
#ifdef USE_DEBUG_LOG
    Log::Indent();
    stringstream s;
    s << "\n####################### Csi2CsiJunction::AddInputCase \n Eliminated vars in order: ";
    for (int i = 0; i < elimOrder.size(); i++)
        s << elimOrder[i].m_var->GetName() << " ";
    s << "\n####################### END AddInputCase\n";
    Log::Write(s.str(), LOG_CSI_JUNCTION);
    Log::Dedent();
#endif

    m_structureUpdated = false;
    m_inMsgsPerCsiCase.push_back(NULL); //inMsgs can be null (in which case no input messages are provided
    std::vector<ConvVar>* v = new std::vector<ConvVar>;
    *v = elimOrder;
    m_elimOrderPerInCase.push_back(v);
    std::vector<ConvVar>* w = new std::vector<ConvVar>;
    m_variablesInInMsgPerInCase.push_back(w); //empty   
}

void Csi2CsiJunction::ComputeWeights()
{
    assert(m_structureUpdated);
    for (int i = 0; i < m_nets.size(); i++)
    {
        for (int j = 0; j < m_nets[i]->size(); j++)
            m_nets[i]->at(j)->ComputeWeights();
    }
    for (int i = 0; i < m_outCsiMsgs.size(); i++)
    {
        if (m_outCsiMsgs[i] != NULL)
            m_outCsiMsgs[i]->InitWeights();
    }
}

void Csi2CsiJunction::Evaluate() //must be fast
{
    assert(m_structureUpdated);
    for (int i = 0; i < m_nets.size(); i++)
    {
        for (int j = 0; j < m_nets[i]->size(); j++)
            m_nets[i]->at(j)->Evaluate();
    }
    for (int i = 0; i < m_prodsPerOutC.size(); i++)
    {
        for (int j = 0; j < m_prodsPerOutC[i]->size(); j++)
        {
            if (m_prodsPerOutC[i] != NULL)
                m_prodsPerOutC[i]->at(j)->Evaluate();
        }
    }
    for (int i = 0; i < m_outCsiMsgs.size(); i++)
    {
        if (m_outCsiMsgs[i] != NULL)
            m_outCsiMsgs[i]->Evaluate();
    }
}

void Csi2CsiJunction::IndicatorsSumOut() //must be fast
{
    assert(m_structureUpdated);
    for (int i = 0; i < m_nets.size(); i++)
    {
        for (int j = 0; j < m_nets[i]->size(); j++)
            m_nets[i]->at(j)->IndicatorsSumOut();
    }
    for (int i = 0; i < m_prodsPerOutC.size(); i++)
    {
        for (int j = 0; j < m_prodsPerOutC[i]->size(); j++)
        {
            if (m_prodsPerOutC[i] != NULL)
                m_prodsPerOutC[i]->at(j)->IndicatorSumOut();
        }
    }
}

void Csi2CsiJunction::IndicatorsAssign(const SPN_VarSet &vars) //must be fast 
{
    assert(m_structureUpdated);
    for (int i = 0; i < m_nets.size(); i++)
    {
        for (int j = 0; j < m_nets[i]->size(); j++)
            m_nets[i]->at(j)->IndicatorsAssign(vars);
    }
    ProdMessage* pm;
    const vector<int>& vals = vars.GetValues();
    for (int k = 0; k < m_prodsPerOutC.size(); k++)
    {
        vector<ProdMessage*>& prods = *(m_prodsPerOutC[k]);
        for (int i = 0; i < prods.size(); i++)
        {
            pm = prods[i];
            int ind = vars.IndexOf(pm->GetIndicators());
            if (ind >= 0)
                pm->IndicatorAssign(vals[ind]); //todo check that it assigns the right one
            else
                pm->IndicatorSumOut();
        }
    }
}

void Csi2CsiJunction::Build()
{
#ifdef USE_DEBUG_LOG
    Log::Indent();
    Log::Write("\n###############################\n############### Csi2CsiJunction::Build() \n", LOG_CSI_JUNCTION);
#endif

    assert(m_structureUpdated == false);
    assert(m_nets.empty());
    assert(GetNInputCases() != 0);
    assert(m_prodsPerOutC.empty());
    assert(!m_csiVarPrefix.empty());

#ifdef USE_DEBUG_LOG
    stringstream msg;
    Log::Write("### THESE ARE THE INPUT MESSAGES: \n", LOG_CSI_JUNCTION);
    for (int in = 0; in < m_inMsgsPerCsiCase.size(); in++)
    {
        stringstream s;
        s << "    Messages for case " << in << ":\n";
        std::vector<SumMessageCommon*>* messages = m_inMsgsPerCsiCase[in];
        if (messages == NULL)
        {
            s << "No messages for this case";
        }
        else
        {
            for (int m = 0; m < messages->size(); m++)
            {
                s << messages->at(m)->ToString() << endl;
            }
        }
        Log::Write(s.str(), LOG_CSI_JUNCTION);
        Log::Write("\n", LOG_CSI_JUNCTION);
    }
    Log::Write("\n", LOG_CSI_JUNCTION);
#endif

    //create a new network from each couple input and output case. 
    //store them in a vector m_nets[i][j] = net with out case i, in case j
    int nOutCases = NOutMessages();
    int nInCases = m_inMsgsPerCsiCase.size();
    for (int outc = 0; outc < nOutCases; outc++)
    {
        vector<SymmetricSPN*>* netsForOutCaseC = new vector<SymmetricSPN*>;

        for (int in = 0; in < nInCases; in++)
        {

            SymmetricSPN* net = new SymmetricSPN;

#ifdef USE_DEBUG_LOG  
            msg.str("");
            msg << "Csi2CsiJunction " << this << " created Symmetric network " << net << " for output case " << outc << " in case " << in;
            Log::Write(msg.str());
#endif

            const vector<FactorFunction*>& factorsForCase = m_factorCreator->GetFactors(*m_variablesInInMsgPerInCase[in], *(m_variablesPerOutCase[outc]));

            net->SetFactors(factorsForCase);
            assert(!m_elimOrderPerInCase[in]->empty());
            net->SetEliminationOrder(ConvertVectorConvVar2SpnVar(*m_elimOrderPerInCase[in]));

            if (m_inMsgsPerCsiCase[in] != NULL)
            {
                const vector<SumMessageCommon*>& inmsgs = *m_inMsgsPerCsiCase[in];
                for (int k = 0; k < inmsgs.size(); k++)
                {
                    assert(!inmsgs[k]->IsRoot() && "the spn root should problably not be in a Csi2CsiJunction. "
                           "Probably one of the input nets was not assigned all the outgoing factors that it needs.");
                    net->AddInput(inmsgs[k]);
                }
            }

#ifdef USE_DEBUG_LOG  
            msg.str("");
            msg << "Csi2CsiJunction " << this << " builds Symmetric network " << net << " for output case " << outc << " in case " << in;
            Log::Write(msg.str());
#endif

            net->Build();
            netsForOutCaseC->push_back(net);
        }
        m_nets.push_back(netsForOutCaseC);

        //now, we create the CSI node prepared for this out case. 
        //for doing this, we create one product message per each symm net case referring to output C
        //and then we pass them to the csi node referring to output C

        vector<ProdMessage*>* prodsForC = new vector<ProdMessage*>;
        for (int in = 0; in < nInCases; in++)
        {
            vector<SumMessageCommon*> inSums;
            const vector<SumMessageCommon*>& outMsgs = netsForOutCaseC->at(in)->GetOutputMessages();
            for (int j = 0; j < outMsgs.size(); j++)
            {
                //                cout<<outMsgs[j]->ToString()<<endl;
                if (outMsgs[j]->IsRoot())
                {
#ifdef USE_DEBUG_LOG  
                    msg.str("");
                    msg << "\nWARNING the network correponding to input case " << in << " and out case " << outc << " has no parents in this case. "
                            "\nWe do not add it as a child. TODO prune it?\n";
                    Log::Write(msg.str());
#endif
                    continue; //todo check. Should we prune the messages that have no parent?
                }
                inSums.push_back(outMsgs[j]);
            }

            if (!inSums.empty())
            {
                ProdMessage* newProd = new ProdMessage(m_nextEliminatedPerOutC[outc], inSums);
                prodsForC->push_back(newProd);
            }
            else
            {
#ifdef USE_DEBUG_LOG  
                msg.str("");
                msg << "\nWARNING no outgoing messages for input case " << in << " and out case " << outc << endl;
                Log::Write(msg.str());
#endif
            }
        }

        if (!prodsForC->empty())
        {

            m_prodsPerOutC.push_back(prodsForC);

            //the csi var of the particular out message has as many states as active outputs.
            stringstream csiVarName;
            csiVarName << m_csiVarPrefix << "_" << outc;
            SPN_Var newCsiVar = SPN_VarManager::PushVar(csiVarName.str(), prodsForC->size());
            //create the factor
            SPN_VarSet csiVarset;
            csiVarset.InsertVariable(newCsiVar);
            TableFactorFunction* newCsiFactor = new TableFactorFunction(csiVarset);
            newCsiFactor->SetUniform();
            m_csiFactors.push_back(newCsiFactor);
            CSISumMessage* csiForOutC = new CSISumMessage(*prodsForC, newCsiFactor, newCsiVar); //todo usage of shared factor

#ifdef USE_DEBUG_LOG  
            msg.str("");
            msg << "The CSI junction " << this << " created CSI message " << csiForOutC << " for output case " << outc;
            Log::Write(msg.str());
#endif
            m_outCsiMsgs.push_back(csiForOutC);
        }
        else
        {
            cout << "\nWARNING no outgoing messages for out case " << outc << endl; //we create however null pointers, as placeholders, in order not to lose the ordering 
            m_prodsPerOutC.push_back(NULL);
            m_outCsiMsgs.push_back(NULL);
        }
    }

    m_structureUpdated = true;

#ifdef USE_DEBUG_LOG
    Log::Write("\n####### END OF Csi2CsiJunction::Build() \n", LOG_CSI_JUNCTION);
    Log::Dedent();
#endif
}
