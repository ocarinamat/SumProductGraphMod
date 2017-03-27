#include "SymmetricSPN.h"

using namespace std;

SymmetricSPN::SymmetricSPN()
: m_structureUpdated(false)
{
}

SymmetricSPN::~SymmetricSPN()
{
    Clean();
}

void SymmetricSPN::Clean()
{
    m_structureUpdated = false;
    for (int i = 0; i < m_prods.size(); i++)
        delete m_prods[i];
    for (int i = 0; i < m_sums.size(); i++)
        delete m_sums[i];
}

void SymmetricSPN::ComputeWeights()
{
    assert(m_structureUpdated);
    int Ns = m_sums.size();
    for (int i = 0; i < Ns; i++)
    {
        m_sums[i]->InitWeights();
    }
}

void SymmetricSPN::SetEliminationOrder(const std::vector<SPN_Var>& varsOrder)
{
    m_eliminationOrder.clear();
    for (int i = 0; i < varsOrder.size(); i++)
    {
        AddEliminationInOrder(varsOrder[i]);
    }
}

void SymmetricSPN::SetFactors(const std::vector<FactorFunction*>& factors)
{
    m_factors.clear();
    for (int i = 0; i < factors.size(); i++)
    {
        AddFactor(*factors[i]);
    }
}

void SymmetricSPN::Evaluate()
{
    assert(m_structureUpdated);

    int Np = m_prods.size();
    int Ns = m_sums.size();
    assert(Np == Ns);

    for (int i = 0; i < Ns; i++)
    {
        m_prods[i]->Evaluate();
        m_sums[i]->Evaluate();
    }
}

void SymmetricSPN::ComputeGradient()
{
    assert(m_structureUpdated);

    int Np = m_prods.size();
    int Ns = m_sums.size();
    assert(Np == Ns);

    for (int i = Ns - 1; i >= 0; i--)
    {
        //        computation in top down order, from root to leaves
        m_sums[i]->ComputeGrad();
        m_prods[i]->ComputeGrad();
    }
}

void SymmetricSPN::IndicatorsAssign(const SPN_VarSet& varset) //this one is done in inference phase, and must be fast
{
    assert(m_structureUpdated);
    ProdMessage* pm;
    const vector<int>& vals = varset.GetValues();

    for (int i = 0; i < m_prods.size(); i++)
    {
        pm = m_prods[i];
        int ind = varset.IndexOf(pm->GetIndicators());
        if (ind >= 0)
            pm->IndicatorAssign(vals[ind]); //todo check that it assigns the right one
        else
            pm->IndicatorSumOut();
    }
}

void SymmetricSPN::IndicatorsSumOut()
{
    IndicatorsAssign(SPN_VarSet());
}

//ProdMessage* SymmetricSPN::GetHead()
//{
//    return m_prods[m_sums.size() - 1];
//}

void SymmetricSPN::AddFactor(FactorFunction& f)
{
    m_structureUpdated = false;

    //check that there are no duplicates
    for (vector<FactorFunction*>::iterator it = m_factors.begin(); it != m_factors.end(); it++)
        assert((*it) != &f);

    m_factors.push_back(&f);
}

void SymmetricSPN::AddEliminationInOrder(const SPN_Var& v)
{
    m_structureUpdated = false;

    //check that it is not already eliminated by any message TODO

    //check that there are no duplicates
    for (int i = 0; i < m_eliminationOrder.size(); i++)
        assert(!v.IsEqual(m_eliminationOrder[i]) && "one of the variables is already in the eliomination order");

    m_eliminationOrder.push_back(v);
}

vector<FactorFunction*> SymmetricSPN::FindAndEliminateFactors(const SPN_Var& v, vector<FactorFunction*>& activeFactors)
{
    vector<FactorFunction*> factorsOnV;
    for (vector<FactorFunction*>::iterator it = activeFactors.begin(); it != activeFactors.end(); it++)
    {
        if ((*it)->GetVarSet().IndexOf(v) >= 0)
        {
            factorsOnV.push_back(*it);
            activeFactors.erase(it);
            it--;
        }
    }
    return factorsOnV;
}

void SymmetricSPN::AddInput(SumMessageCommon* sm)
{
    for (int i = 0; i < m_unusedSums.size(); i++)
    {
        assert(m_unusedSums[i] != sm);
    }
    m_unusedSums.push_back(sm);
}

int SymmetricSPN::ComputeMaxWidth()
{
    assert(m_structureUpdated);
    int Np = m_prods.size();
    assert(Np == m_sums.size());

    int max = -1;
    int curr = -1;
    for (int i = 0; i < Np; i++)
    {
        curr = m_prods[i]->ComputeWidth();
        if (curr > max)
            max = curr;
        curr = m_sums[i]->ComputeWidth();
        if (curr > max)
            max = curr;
    }
    return max;
}

std::vector<SumMessageCommon*> SymmetricSPN::FindAndEliminateMsgByVariables(const SPN_Var& v)
{
    vector<SumMessageCommon*> sumsOnV;
    for (vector<SumMessageCommon*>::iterator it = m_unusedSums.begin(); it != m_unusedSums.end(); it++)
    {
        if ((*it)->GetOutVars().IndexOf(v) >= 0) //todo efficiency
        {
            sumsOnV.push_back(*it);
            m_unusedSums.erase(it);
            it--;
        }
    }
    return sumsOnV;
}

const SumMessageCommon* SymmetricSPN::GetRoot() const
{
    assert(m_unusedSums.size() == 1);
    assert(m_structureUpdated);
    return m_unusedSums[0];
}

string SymmetricSPN::ToString() const
{
    stringstream s;
    s << "###########  SymmetricSPN INFO: ID " << this << "\n";
    s << "Input Messages:\n";
    for (int f = 0; f < m_unusedSums.size(); f++)
        s << "    " << m_unusedSums[f]->ToString() << "\n";
    s << "Factors:\n";
    for (int f = 0; f < m_factors.size(); f++)
        s << "    " << m_factors[f]->ToStringCompact() << "\n";
    s << "\nelimination order:\n";
    for (int f = 0; f < m_eliminationOrder.size(); f++)
        s << "    " << m_eliminationOrder[f].GetName() << ", ";
    s << "\n\nEND OF SymmetricSPN INFO #######";
    return s.str();
}

void SymmetricSPN::Build()
{
#ifdef USE_DEBUG_LOG
    Log::Indent();
    Log::Write("\n#######################", LOG_SYMM_SPN);
    Log::Write("SymmetricSPN::Build()", LOG_SYMM_SPN);
    Log::Write(ToString(), LOG_SYMM_SPN);
    Log::Write("\n\n", LOG_SYMM_SPN);
    stringstream msg;
#endif

    assert(!m_factors.empty());
    //take a copy because we modify it
    vector<FactorFunction*> activeFactors = m_factors;

    Clean();

    ProdMessage* prod = NULL;
    SumMessage* sum = NULL;
    int Nvars = m_eliminationOrder.size();

    for (int i = 0; i < Nvars; i++)
    {
        const SPN_Var& v = m_eliminationOrder[i];

#ifdef USE_DEBUG_LOG
        Log::Indent();
        Log::Write("ELIMINATING ", LOG_SYMM_SPN_DETAILS);
        Log::Write(v.GetName(), LOG_SYMM_SPN_DETAILS);
        Log::Endl(LOG_SYMM_SPN_DETAILS);
#endif

        /////////////////////////////////////////////
        //we create a product message in which the indicator variable (the one currently eliminated) is v itsel
        //and it takes as input the product of all the messages that depend on that
        /////////////////////////////////////////////

        const std::vector<SumMessageCommon*> &inSums = FindAndEliminateMsgByVariables(v);

        if (inSums.empty())
        {
            prod = new ProdMessage(v);
        }
        else
        {
            prod = new ProdMessage(v, inSums);
        }

        //add to the array of prod messages
        m_prods.push_back(prod);


#ifdef USE_DEBUG_LOG
        msg.str("");
        msg << "created product msg " << prod << " : " << prod->ToString();
        Log::Write(msg.str());
        Log::Write("that multiplies: ", LOG_SYMM_SPN_DETAILS);
        for (int ll = 0; ll < inSums.size(); ll++)
        {
            Log::Write("\t", LOG_SYMM_SPN_DETAILS);
            Log::Write(inSums[ll]->ToString(), LOG_SYMM_SPN_DETAILS);
        }
        Log::Endl(LOG_SYMM_SPN_DETAILS);
#endif

        /////////////////////////////////////////////
        //we then create a sum message. 
        //the factors that we use is the product of all the ones in which there is v
        /////////////////////////////////////////////
        //select the factors...
        const vector<FactorFunction*>& selectedFactors = FindAndEliminateFactors(v,activeFactors);

        //create the sum message
        sum = new SumMessage(*prod, selectedFactors);

        //TODO? the output variables are the variables that appear in any factor with v, plus the ones that are already there    
        m_unusedSums.push_back(sum);
        m_sums.push_back(sum);

#ifdef USE_DEBUG_LOG
        msg.str("");
        msg << "created sum msg " << sum << " : " << sum->ToString();
        Log::Write(msg.str());
        Log::Endl(LOG_SYMM_SPN_DETAILS);
#endif
#ifdef USE_DEBUG_LOG
        Log::Write("Using these factors: ", LOG_SYMM_SPN_DETAILS);
        for (int ff = 0; ff < selectedFactors.size(); ff++)
            Log::Write(selectedFactors[ff]->ToStringCompact(), LOG_SYMM_SPN_DETAILS);
        Log::Endl(LOG_SYMM_SPN_DETAILS);
        Log::Dedent();
#endif


    }
    m_structureUpdated = true;

#ifdef USE_DEBUG_LOG
    Log::Write("CREATED THESE OUTPUT MESSAGES: ", LOG_SYMM_SPN);
    Log::Indent();
    for (int i = 0; i < GetOutputMessages().size(); i++)
        Log::Write(GetOutputMessages()[i]->ToString(), LOG_SYMM_SPN);
    Log::Dedent();

    Log::Write("#### end of SYmmSPN::Build ###################", LOG_SYMM_SPN);
    Log::Dedent();
#endif
}

bool SymmetricSPN::Eliminates(const SPN_Var& v) const
{
    for (int i = 0; i < m_eliminationOrder.size(); i++)
    {
        if (v.IsEqual(m_eliminationOrder[i]))
            return true;
    }
    return false;
}

void SymmetricSPN::CheckGradient_Weights()
{
    cout << "\n\nCheckGradient_Weights\n\n";
    vector<double> SampleGradient;
    vector<double> AnalyticalGradient;
    double eps = 0.000001;

    //first compute the analytical gradient
    ComputeGradient();

    //select every sum node in every sum message
    for (int s = m_sums.size() - 1; s >= 0; s--)
    {
        //select every node
        int nNodes = m_sums[s]->ComputeWidth();
        for (int n = 0; n < nNodes; n++)
        {
            SumNode& currNode = m_sums[s]->GetNode(n);
            //evaluate the model changing one weight at a time
            const double* originalWeights = currNode.GetWeights();
            vector<double> weights(originalWeights, originalWeights + currNode.GetNChildren());
            cout << "level " << s << " node " << n << endl;

            for (int w = 0; w < currNode.GetNChildren(); w++)
            {
                double grad = currNode.GetGradient(w);
                double initialVal = originalWeights[w];
                weights[w] = initialVal + eps;
                currNode.SetWeights(weights);
                Evaluate();
                double valPlus = GetRootVal();

                weights[w] = initialVal - eps;
                currNode.SetWeights(weights);
                Evaluate();
                double valMinus = GetRootVal();

                //add the computed values
                double sampleGrad = (valPlus - valMinus) / (2 * eps);
                cout << "    w " << w << ": grad " << grad << " sampleGrad " << sampleGrad << " diff " << grad - sampleGrad << endl;

                //reset the weights
                weights[w] = initialVal;
                currNode.SetWeights(weights);

                SampleGradient.push_back(sampleGrad);
                AnalyticalGradient.push_back(grad);
            }
        }
    }
}

//void SymmetricSPN::CheckGradient_Activation()
//{
//    vector<double> SampleGradient;
//    vector<double> AnalyticalGradient;
//    double eps = 0.000001;
//
//    //first compute the analytical gradient
//    ComputeGradient();
//
//    //select every sum node in every sum message
//    for (int s = 0; s < m_sums.size(); s++)
//    {
//        //select every node
//        int nNodes = m_sums[s]->ComputeWidth();
//        for (int n = 0; n < nNodes; n++)
//        {
//            SumNode& currNode = m_sums[s]->GetNode(n);
//
//            //evaluate the model changing one value at a time
//            double grad = currNode.GetDspnDval();
//            double initialVal = currNode.GetVal();
//
//            currNode.SetVal(initialVal + eps);
//            Evaluate();
//            double valPlus = GetRoot()->GetNode(0).GetVal();
//
//            currNode.SetVal(initialVal - eps);
//            Evaluate();
//            double valMinus = GetRoot()->GetNode(0).GetVal();
//
//            //add the computed gradients
//            double sampleGrad = (valPlus - valMinus) / (2 * eps);
//            cout << "grad " << grad << " sampleGrad " << sampleGrad << " diff " << grad - sampleGrad << endl;
//
//            //reset the value
//            currNode.SetVal(initialVal);
//
//            SampleGradient.push_back(sampleGrad);
//            AnalyticalGradient.push_back(grad);
//        }
//    }
//
//    //select every node in every product message
//    for (int s = 0; s < m_prods.size(); s++)
//    {
//        //select every node
//        int nNodes = m_prods[s]->ComputeWidth();
//        for (int n = 0; n < nNodes; n++)
//        {
//            ProdNode& currNode = m_prods[s]->GetNode(n);
//
//            //evaluate the model changing one value at a time
//            double grad = currNode.GetDspnDval();
//            double initialVal = currNode.GetVal();
//
//            currNode.SetVal(initialVal + eps);
//            Evaluate();
//            double valPlus = GetRoot()->GetNode(0).GetVal();
//
//            currNode.SetVal(initialVal - eps);
//            Evaluate();
//            double valMinus = GetRoot()->GetNode(0).GetVal();
//
//            //add the computed gradients
//            double sampleGrad = (valPlus - valMinus) / (2 * eps);
//            cout << "grad " << grad << " sampleGrad " << sampleGrad << " diff " << grad - sampleGrad << endl;
//
//            //reset the value
//            currNode.SetVal(initialVal);
//
//            SampleGradient.push_back(sampleGrad);
//            AnalyticalGradient.push_back(grad);
//        }
//    }
//}