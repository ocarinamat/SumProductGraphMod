#include "SPGMnodes.h"
#include <iostream>
#include <map>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace std;

SumNode::SumNode(Vertex v)
: SPGMnode(SPGMnode::sum, v)
{
}

ProdNode::ProdNode(Vertex v) : SPGMnode(SPGMnode::prod, v)
{
}

std::string SumNode::ToStringWeights()
{
    std::stringstream str;
    str << this->ToString() << "weights: [";
    for (int i = 0; i < m_weights.size(); i++)
    {
        str << m_weights[i] << ", ";
    }
    str << "]";
    return str.str();
}

std::string SPGMnode::ToStringMessages(bool printParams)
{
    std::stringstream str;
    str << "Messages of node " << ToString() << ":\n";
    for (VertexMsgMap::iterator it = m_messages.begin(); it != m_messages.end(); it++)
    {
        Message* msg = *it;
        Vertex toVertex = (msg->m_toVertex);
        str << "  msg " <<msg<<" from " << m_vertex << " to " << toVertex;
        if (m_type == vnode)
        {
            if (printParams)
                str << "  with factor:\n " << msg->m_factor->ToString();
            else
                str << "  with factor " << msg->m_factor;
        }
        str << "\n ch messages: ";
        for (int i=0; i<msg->m_childMsgs.size(); i++)
        {
            str <<msg->m_childMsgs[i]<<", ";
        }
        str << "\n  val: " << msg->m_logMessage.ToStringExponentiated();
        str << std::endl;
    }
    return str.str();
}

void SumNode::SetWeights(const std::vector<Real>& w)
{
    m_weights = w;

    Real sum = 0;
    for (int i = 0; i < w.size(); i++)
    {
        sum += w[i];
        assert(!isnan(w[i]));
        assert(w[i] >= 0);
    }

    if (abs(sum - 1.f) > 0.0000000000001)
    {
        std::stringstream errm;
        errm << "\nError:\nassigned weights are not normalized\n";
        errm << "sum is " << sum << " and the weights are: \n";

        for (int i = 0; i < w.size(); i++)
            errm << w[i] << ", ";
        std::cerr << errm.str() << std::endl;
    }

    for (int i = 0; i < m_weights.size(); i++)
    {
        m_weights[i] /= sum;
    }

    ResetBeta();
}

void SPGMnode::CreateMessages(int N, int dim)
{
    for (VertexMsgMap::iterator it = m_messages.begin(); it != m_messages.end(); it++)
    {
        Message* msg = *it;
        if (msg->m_toVertex == ROOT_VERTEX)
        {
            msg->m_logMessage = RealMatrix(N, 1);
            msg->m_logDSpnDmessage = RealMatrix(N, 1);
        }
        else
        {
            msg->m_logMessage = RealMatrix(N, dim);
            msg->m_logDSpnDmessage = RealMatrix(N, dim);
        }
    }
}

void SumNode::SetWeightsRandom()
{
    m_weights = std::vector<Real>(m_messages[0]->m_childMsgs.size(), std::numeric_limits<Real>::quiet_NaN());

    Real sum = 0;
    for (int i = 0; i < m_weights.size(); i++)
    {
        Real rv = SPN_Rand(RAND_INIT_WEIGHT_MIN, RAND_INIT_WEIGHT_MAX);
        sum += rv;
        m_weights[i] = rv;
    }

    Real invSum = 1.f / sum;

    for (int i = 0; i < m_weights.size(); i++)
        m_weights[i] *= invSum;
}

void Vnode::Eval(const DataMatrix& data)
{
    bool verbose = false;
    assert(data.nVars() > m_var);
    int N = data.nData();

    for (VertexMsgMap::iterator it = m_messages.begin(); it != m_messages.end(); it++)
    {
        Message* msg = *it;
        RealMatrix& logMessage = (msg->m_logMessage);
        CPTfactor* factor = msg->m_factor;
        int toVar = msg->m_toVar;
        Vertex toVertex = msg->m_toVertex;

        //reset derivative
        RealMatrix& logDSpnDmessage = (msg->m_logDSpnDmessage);
        logDSpnDmessage.SetVal(ZERO_LOGVAL);

        //eval all the outgoing messages
        vector<Message*>& currChildMessages = msg->m_childMsgs;
        assert(currChildMessages.size() <= 1);


        if (verbose)
            std::printf("evaluating message from %d to %d with currState %s parentState %s\n", (int) m_vertex, (int) toVertex, VecToString<byte>(data.ByCol().GetRow(m_var)).c_str(), VecToString<byte>(data.ByCol().GetRow(toVar)).c_str());

        if (!currChildMessages.empty())
            factor->Eval((currChildMessages[0]->m_logMessage), logMessage, data, m_var, toVar);
        else
            factor->Eval(logMessage, data, m_var, toVar);

    }
}

void SumNode::Eval()
{
    int N = m_messages.front()->m_logMessage.nRows();
    int dim = m_messages.front()->m_logMessage.nCols();
    int nChildren = m_messages.front()->m_childMsgs.size();
    RealMatrix maxChVal;

    std::vector<Real> logW(m_weights.size());
    for (int i = 0; i < m_weights.size(); i++)
        logW[i] = log(m_weights[i]);

    //    cout<<VecToStringExponentiated(logW)<<endl;

    //for each output message... 
    for (VertexMsgMap::iterator it = m_messages.begin(); it != m_messages.end(); it++)
    {
        /*evaluates the element by element weighted sum of the input vectors, which are under log:
         * if input vectors log A, log B, it computes log(w1 A + w2 B)
         */
        Message* msg = *it;
        vector<Message*>& currChildMessages = msg->m_childMsgs;
        assert(m_weights.size() == currChildMessages.size());
        RealMatrix& logMessage = (msg->m_logMessage);
        RealMatrix& logDSpnDmessage = (msg->m_logDSpnDmessage);
        Vertex toVertex = msg->m_toVertex;

        logDSpnDmessage.SetVal(ZERO_LOGVAL);

        //find max ch val for each entry, between all possible children, 
        maxChVal = (currChildMessages[0]->m_logMessage);
        for (int i = 1; i < nChildren; i++) // start from second child
        {
            RealMatrix& chMsg = (currChildMessages[i]->m_logMessage);
            for (int n = 0; n < N; n++)
            {
                for (int j = 0; j < dim; j++)
                {
                    Real chVal = chMsg.GetVal(n, j);
                    Real& mv = maxChVal.EditVal(n, j);
                    if (chVal > mv)
                    {
                        mv = chVal;
                        //#ifndef NDEBUG
                        //                        if (isnan(mv))
                        //                        {
                        //                            std::printf("max val is nan");
                        //                            assert(false);
                        //                        }
                        //#endif
                    }
                }
            }
        }

        //        cout<<"maxChVal : \n"<<maxChVal.ToStringExponentiated()<<endl;

        //now compute the new message
        logMessage.SetVal(0);
        for (int i = 0; i < nChildren; i++)
        {
            Real logWi = logW[i];
            const RealMatrix& chMsg = (currChildMessages[i]->m_logMessage);

            //                    cout<<"child msg: \n"<<chMsg.ToStringExponentiated()<<endl;

            for (int n = 0; n < N; n++)
            {
                for (int j = 0; j < dim; j++)
                {
                    Real maxChV = maxChVal.GetVal(n, j);
                    if (!isinf(maxChV))
                    {
                        Real tempVal = exp(chMsg.GetVal(n, j) - maxChVal.GetVal(n, j) + logWi);
                        Real& logJN = logMessage.EditVal(n, j);
                        logJN += tempVal;
#ifndef NDEBUG
                        if (isnan(tempVal))
                        {
                            std::printf("logMessage has nan");
                            assert(false);
                        }
#endif
                    }
                }
            }
        }
        //        logMessage.ApplyLog();
        //        logMessage.Sum(maxChVal);
        for (int n = 0; n < N; n++)
        {
            for (int j = 0; j < dim; j++)
            {
                Real maxChV = maxChVal.GetVal(n, j);
                Real& logJN = logMessage.EditVal(n, j);
                if (isinf(maxChV))
                {
                    logJN = maxChV; //this avoids problems when maxChV == -inf
                }
                else
                {
                    logJN = log(logJN) + maxChV;
#ifndef NDEBUG
                    if (isnan(logJN))
                    {
                        std::printf("logMessage has nan");
                        assert(false);
                    }
#endif
                }
            }
        }
        //        cout<<"produced Message: \n"<<logMessage.ToStringExponentiated()<<endl;
    }

}

void SumNode::StackBeta(const vector<Real>& sampleLogProb)
{
    int N = sampleLogProb.size();
    for (VertexMsgMap::iterator it = m_messages.begin(); it != m_messages.end(); it++)
    {
        /*evaluates the element by element weighted sum of the input vectors, which are under log:
         * if input vectors log A, log B, it computes log(w1 A + w2 B)
         */
        Message* msg = *it;
        vector<Message*>& currChildMessages = msg->m_childMsgs;
        //        std::vector<Real>& logMessage = *(msg->m_logMessage);
        const RealMatrix& logDSpnDmessage = (msg->m_logDSpnDmessage);
        Vertex toVertex = msg->m_toVertex;

        for (int i = 0; i < m_weights.size(); i++)
        {
            assert(currChildMessages[i]->m_toVertex == toVertex);
            const RealMatrix& chMsg = currChildMessages[i]->m_logMessage;
            for (int n = 0; n < N; n++)
            {
                Real sampleLL = sampleLogProb[n];
                for (int j = 0; j < chMsg.nRows(); j++)
                {
                    Real logBetaQJN = logDSpnDmessage.GetVal(n, j) + chMsg.GetVal(n, j) - sampleLL;
                    if (logBetaQJN != ZERO_LOGVAL)
                        m_logBeta[i] = AddLog(m_logBeta[i], logBetaQJN);
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

void SumNode::PassDerivative_andStackBeta(const vector<Real>& sampleLogProb)
{
    int N = m_messages[0]->m_logMessage.nRows();
    for (VertexMsgMap::iterator it = m_messages.begin(); it != m_messages.end(); it++)
    {
        /*evaluates the element by element weighted sum of the input vectors, which are under log:
         * if input vectors log A, log B, it computes log(w1 A + w2 B)
         */
        Message* msg = *it;
        vector<Message*>& currChildMessages = msg->m_childMsgs;
        RealMatrix& logDSpnDmessage = (msg->m_logDSpnDmessage);

        for (int i = 0; i < m_weights.size(); i++)
        {
            assert(currChildMessages[i]->m_toVertex == msg->m_toVertex);

            Real w = m_weights[i];
            RealMatrix& chDer = (currChildMessages[i]->m_logDSpnDmessage);
            const RealMatrix& chMsg = currChildMessages[i]->m_logMessage;

            for (int n = 0; n < N; n++)
            {
                for (int j = 0; j < logDSpnDmessage.nCols(); j++)
                {
                    Real sampleLL = sampleLogProb[n];

                    Real logDerivative = logDSpnDmessage.GetVal(n, j);

                    //                    assert(!isnan(chDer.GetVal(n,j)));

                    if (logDerivative != ZERO_LOGVAL)
                    {
                        //update input derivative
                        Real l = logDerivative + log(w);
                        Real& chDerVal = chDer.EditVal(n, j);
                        if (chDerVal == ZERO_LOGVAL)
                            chDerVal = l;
                        else
                            chDerVal = AddLog(chDerVal, l);
                        //////////////////////////
                        // stack beta
                        Real logBetaQJN = logDerivative + chMsg.GetVal(n, j) - sampleLL;
                        if (logBetaQJN != ZERO_LOGVAL)
                            m_logBeta[i] = AddLog(m_logBeta[i], logBetaQJN);
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
}

void SumNode::PassDerivative()
{
    int N = m_messages[0]->m_logMessage.nRows();
    for (VertexMsgMap::iterator it = m_messages.begin(); it != m_messages.end(); it++)
    {
        /*evaluates the element by element weighted sum of the input vectors, which are under log:
         * if input vectors log A, log B, it computes log(w1 A + w2 B)
         */
        Message* msg = *it;
        vector<Message*>& currChildMessages = msg->m_childMsgs;
        RealMatrix& logDSpnDmessage = (msg->m_logDSpnDmessage);

        for (int i = 0; i < m_weights.size(); i++)
        {
            assert(currChildMessages[i]->m_toVertex == msg->m_toVertex);

            Real w = m_weights[i];
            RealMatrix& chDer = (currChildMessages[i]->m_logDSpnDmessage);

            for (int n = 0; n < N; n++)
            {
                for (int j = 0; j < logDSpnDmessage.nCols(); j++)
                {
                    Real logDerivative = logDSpnDmessage.GetVal(n, j);

                    //                    assert(!isnan(chDer.GetVal(n,j)));

                    if (logDerivative != ZERO_LOGVAL)
                    {
                        //update input derivative
                        Real l = logDerivative + log(w);
                        Real& chDerVal = chDer.EditVal(n, j);
                        if (chDerVal == ZERO_LOGVAL)
                            chDerVal = l;
                        else
                            chDerVal = AddLog(chDerVal, l);
                        //////////////////////////
                    }
                }
            }
        }
    }
}

void ProdNode::PassDerivative()
{
    int N = m_messages[0]->m_logMessage.nRows();
    for (VertexMsgMap::iterator it = m_messages.begin(); it != m_messages.end(); it++)
    {
        /*evaluates the element by element weighted sum of the input vectors, which are under log:
         * if input vectors log A, log B, it computes log(w1 A + w2 B)
         */
        Message* msg = *it;
        vector<Message*>& currChildMessages = msg->m_childMsgs;
        RealMatrix& m_logMessage = (msg->m_logMessage);
        RealMatrix& m_logDSpnDmessage = (msg->m_logDSpnDmessage);

        for (int i = 0; i < currChildMessages.size(); i++)
        {
            assert(currChildMessages[i]->m_toVertex == msg->m_toVertex);

            RealMatrix& chDer = (currChildMessages[i]->m_logDSpnDmessage);
            RealMatrix& chMsg = (currChildMessages[i]->m_logMessage);

            for (int n = 0; n < N; n++)
            {
                for (int j = 0; j < m_logMessage.nCols(); j++)
                {
                    Real logDerivative = m_logDSpnDmessage.GetVal(n, j);

                    if (logDerivative != ZERO_LOGVAL)
                    {
                        //compute current input derivative
                        Real l;
                        Real chVal = chMsg.GetVal(n, j);
                        if (chVal == ZERO_LOGVAL) //the derivative is the product (log sum) of all other elements 
                        {
                            l = 0;
                            for (int otherCh = 0; otherCh < currChildMessages.size(); otherCh++)
                                if (otherCh != j)
                                    l += chMsg.GetVal(n, otherCh);
                        }
                        else
                        {
                            l = logDerivative + m_logMessage.GetVal(n, j) - chVal;
                        }

                        //add it 
                        Real& chDer_curr = chDer.EditVal(n, j);
                        if (chDer_curr == ZERO_LOGVAL)
                            chDer_curr = l;
                        else
                            chDer_curr = AddLog(chDer_curr, l);
                    }
                    //////////////////////////
                }
            }
        }
    }

}

//void ProdNode::PassDerivative()
//{
//    int N = m_messages[0]->m_logMessage.nCols();
//    for (VertexMsgMap::iterator it = m_messages.begin(); it != m_messages.end(); it++)
//    {
//        /*evaluates the element by element weighted sum of the input vectors, which are under log:
//         * if input vectors log A, log B, it computes log(w1 A + w2 B)
//         */
//        Message* msg = *it;
//        vector<Message*>& currChildMessages = msg->m_childMsgs;
//        RealMatrix& m_logMessage = (msg->m_logMessage);
//        RealMatrix& m_logDSpnDmessage = (msg->m_logDSpnDmessage);
//
//        for (int i = 0; i < currChildMessages.size(); i++)
//        {
//            assert(currChildMessages[i]->m_toVertex == msg->m_toVertex);
//
//            RealMatrix& chDer = (currChildMessages[i]->m_logDSpnDmessage);
//            RealMatrix& chMsg = (currChildMessages[i]->m_logMessage);
//
//            for (int j = 0; j < m_logMessage.nRows(); j++)
//            {
//                for (int n = 0; n < N; n++)
//                {
//                    Real logDerivative = m_logDSpnDmessage.GetVal(n,j);
//
//                    if (logDerivative != ZERO_LOGVAL)
//                    {
//                        //compute current input derivative
//                        Real l;
//                        Real chVal = chMsg.GetVal(n,j);
//                        if (chVal == ZERO_LOGVAL) //the derivative is the product (log sum) of all other elements 
//                        {
//                            l = 0;
//                            for (int otherCh = 0; otherCh < currChildMessages.size(); otherCh++)
//                                if (otherCh != j)
//                                    l += chMsg.GetVal(otherCh, n);
//                        }
//                        else
//                        {
//                            l = logDerivative + m_logMessage.GetVal(n,j) - chMsg.GetVal(n,j);
//                        }
//
//                        //add it 
//                        Real& chDer_curr = chDer.EditVal(n,j);
//                        if (chDer_curr == ZERO_LOGVAL)
//                            chDer_curr = l;
//                        else
//                            chDer_curr = AddLog(chDer_curr, l);
//                    }
//                    //////////////////////////
//                }
//            }
//        }
//    }
//
//}

void ProdNode::Eval()
{
    for (VertexMsgMap::iterator it = m_messages.begin(); it != m_messages.end(); it++)
    {
        /*evaluates the element by element weighted sum of the input vectors, which are under log:
         * if input vectors log A, log B, it computes log(w1 A + w2 B)
         */
        Message* msg = *it;
        vector<Message*>& currChildMessages = msg->m_childMsgs;
        RealMatrix& m_logMessage = msg->m_logMessage;
        RealMatrix& m_logDSpnDmessage = msg->m_logDSpnDmessage;
        int toVertex = msg->m_toVertex;

        /*evaluates the element by element product of the input vectors, which are under log:
         * if input vectors log A, log B, it computes log(A .* B)
         */

        m_logDSpnDmessage.SetVal(ZERO_LOGVAL);
        m_logMessage = (currChildMessages[0]->m_logMessage);

        for (int i = 1; i < currChildMessages.size(); i++)
        {
            assert(currChildMessages[i]->m_toVertex == toVertex);
            RealMatrix& chMsg = (currChildMessages[i]->m_logMessage);
            m_logMessage.Sum(chMsg);
        }

#ifndef NDEBUG
        if (m_logMessage.HasNan())
        {
            std::printf("m_logMessage has nan");
            std::exit(1);
        }
#endif
    }
}

void SumNode::EMupdate(Real uniformPrior)
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

    ResetBeta();

    if (sumTempBeta == 0)
    {
        //no responsibility for this node. return. 
        return;
    }

    assert(uniformPrior>=0);
    int nw = newW.size();
    Real z = sumTempBeta + uniformPrior;
    Real uOverNw = uniformPrior / nw;
    for (int k = 0; k < nw; k++)
    {
        m_weights[k] = (newW[k] + uOverNw)/z;
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

//void Vnode::PassDerivative(const DataMatrix& data)
//{
//    const vector<byte>& dataRow = data.GetRow(m_var);
//    bool verbose = false;
//    if (verbose)
//    {
//        for (VertexMsgMap::iterator it = m_messages.begin(); it != m_messages.end(); it++)
//        {
//            Message* msg = *it;
//            RealMatrix& logDSpnDmessage = (msg->m_logDSpnDmessage);
//            if (verbose)
//                std::printf("passing derivative of message from %d to %d which is %s \n", (int) m_vertex, (int) msg->m_toVertex, logDSpnDmessage.ToStringExponentiated().c_str());
//        }
//    }
//
//    vector<Message*>& currChildMessages = m_messages[0]->m_childMsgs; //it is the SAME for every output message (since the vnode has only one child))
//    assert(currChildMessages.size() <= 1);
//    if (currChildMessages.empty())
//        return;
//
//    for (VertexMsgMap::iterator it = m_messages.begin(); it != m_messages.end(); it++)
//    {
//        Message* msg = *it;
//        RealMatrix& logDSpnDmessage = (msg->m_logDSpnDmessage);
//        CPTfactor* factor = msg->m_factor;
//
//        factor->PassDerivative_andStackBeta(logDSpnDmessage, (currChildMessages[0]->m_logDSpnDmessage), dataRow);
//    }
//}

void Vnode::PassDerivative_andStackBeta(const DataMatrix& data, const std::vector<Real>& sampleLogProb)
{
    const vector<byte>& dataRow = data.ByCol().GetRow(m_var);
    bool verbose = false;
    if (verbose)
    {
        for (VertexMsgMap::iterator it = m_messages.begin(); it != m_messages.end(); it++)
        {
            Message* msg = *it;
            RealMatrix& logDSpnDmessage = (msg->m_logDSpnDmessage);
            if (verbose)
                std::printf("passing derivative of message from %d to %d which is %s \n", (int) m_vertex, (int) msg->m_toVertex, logDSpnDmessage.ToStringExponentiated().c_str());
        }
    }

    vector<Message*>& currChildMessages = m_messages[0]->m_childMsgs; //it is the SAME for every output message (since the vnode has only one child))
    assert(currChildMessages.size() <= 1);
    if (!currChildMessages.empty())
        for (VertexMsgMap::iterator it = m_messages.begin(); it != m_messages.end(); it++)
        {
            Message* msg = *it;
            RealMatrix& logDSpnDmessage = (msg->m_logDSpnDmessage);
            CPTfactor* factor = msg->m_factor;

            factor->PassDerivative_andStackBeta(logDSpnDmessage, (currChildMessages[0]->m_logDSpnDmessage), dataRow, (currChildMessages[0]->m_logMessage), sampleLogProb);
        }
    else
    {
        RealMatrix emptyMat;
        for (VertexMsgMap::iterator it = m_messages.begin(); it != m_messages.end(); it++)
        {
            Message* msg = *it;
            RealMatrix& logDSpnDmessage = (msg->m_logDSpnDmessage);
            CPTfactor* factor = msg->m_factor;

            factor->StackBeta(emptyMat, logDSpnDmessage, dataRow, sampleLogProb);
        }
    }
}