/* 
 * File:   Node.h
 * Author: mdesana
 *
 * Created on October 2016
 */

#pragma once

#include "Common.h"
#include "Factor.h"

class Message
{
public:
    RealMatrix m_logMessage;
    RealMatrix m_logDSpnDmessage;
    CPTfactor* m_factor;
    int m_toVar;
    int m_toVertex;
    std::vector<Message*> m_childMsgs;

    Message() : m_toVar(VAR_NON_INITIALIZED), m_factor(NULL), m_toVertex(EMPTY_VERTEX)
    {
        m_childMsgs = std::vector<Message*>();
    }
    
//    ~Message() 
//    {
//        std::cout<<"deleting message "<<this<<std::endl;
//    }
    
    void Init(int toVar, int toVertex)
    {
        assert(m_toVar == VAR_NON_INITIALIZED);
        assert(m_toVertex == EMPTY_VERTEX);
        m_toVertex = toVertex;
        m_toVar = toVar;
    }
};

typedef std::vector<Message*> VertexMsgMap; //todo use Vertex, not var. 

class SPGMnode
{
public:

    enum type
    {
        prod, sum, vnode
    };

    SPGMnode(type t, Vertex v) : m_type(t), m_vertex(v)
    {
    }

    ~SPGMnode()
    {
//        if (m_verbose)
//            std::cout<<"deleting SPGMnode "<<this<<std::endl;
        ClearMessages();
    }

    void ClearChildren()
    {
        for (int i = 0; i < m_messages.size(); i++)
            m_messages[i]->m_childMsgs.clear();
    }

    void ClearMessages()
    {
        for (VertexMsgMap::iterator it = m_messages.begin(); it != m_messages.end(); it++)
            delete (*it);

        m_messages = VertexMsgMap();
    }

    int NChildren()
    {
        if (m_messages[0]->m_childMsgs.empty())
            return 0;
        else
            return m_messages[0]->m_childMsgs.size(); //the number of children is the same for all messages
    }

    Vertex GetVertex() const
    {
        return m_vertex;
    }

    void CreateMessages(int N, int dim); //initializes messages of the proper size

    void AddChildMsg(Message* m, Vertex toVertex)
    {
        assert(m != NULL);
        Message* mTo = GetMessage(toVertex);
        for (int i = 0; i < mTo->m_childMsgs.size(); i++)
            assert(mTo->m_childMsgs[i] != m && "child message already present");
        mTo->m_childMsgs.push_back(m);
    }

    void CopyData(SPGMnode* other, const std::map<const CPTfactor*, CPTfactor*>& correspondingFactors) const 
    {
        assert(other != NULL);
        other->m_type = m_type;
        other->m_vertex = m_vertex;
        other->m_messages = VertexMsgMap(m_messages.size());
        for (int i = 0; i < m_messages.size(); i++)
        {
            other->m_messages[i] = new Message();
            other->m_messages[i]->Init(m_messages[i]->m_toVar, m_messages[i]->m_toVertex);
            const CPTfactor* thisFactor = m_messages[i]->m_factor;
            CPTfactor* otherFactor = correspondingFactors.at(thisFactor);
            assert(otherFactor!=NULL);
            other->m_messages[i]->m_factor = otherFactor;
        }
    }

    void CopyData(SPGMnode* other) const 
    {
        assert(other != NULL);
        other->m_type = m_type;
        other->m_vertex = m_vertex;
        other->m_messages = VertexMsgMap(m_messages.size());
        for (int i = 0; i < m_messages.size(); i++)
        {
            other->m_messages[i] = new Message();
            other->m_messages[i]->Init(m_messages[i]->m_toVar, m_messages[i]->m_toVertex);
            other->m_messages[i]->m_factor = NULL;
        }
    }

    std::list<Vertex> GetOutputMessageVertices()
    {
        std::list<Vertex> outList;
        for (VertexMsgMap::iterator it = m_messages.begin(); it != m_messages.end(); it++)
        {
            outList.push_back((*it)->m_toVertex);
        }
        return outList;
    }

    const type GetType() const
    {
        return m_type;
    }

    virtual void Eval()
    {
    };

    virtual void PassDerivative() = 0;

    int Nmessages() const
    {
        return m_messages.size();
    }

    RealMatrix& EditMsgLogVals(Vertex toVertex)
    {
        Message* m = GetMessage(toVertex);
        assert(m != NULL);
        return (m->m_logMessage);
    }

    RealMatrix& EditLogDspnDmessage(Vertex toVertex)
    {
        Message* m = GetMessage(toVertex);
        assert(m != NULL);
        return (m->m_logDSpnDmessage);
    }

    Message* GetMessage(Vertex toVertex)
    {
        for (VertexMsgMap::iterator it = m_messages.begin(); it != m_messages.end(); it++)
            if ((*it)->m_toVertex == toVertex)
                return *it;
        return NULL;
    }

    Message& AddMessage(Vertex toVertex, int toVar)
    {
        Message* found = GetMessage(toVertex);
        assert(found == NULL); //message must not be already in the map
        Message* m = new Message();
        m_messages.push_back(m);
        m->Init(toVar, toVertex);
        return *m;
    }

    virtual std::string ToString() = 0;

    std::string ToStringMessages(bool printParams);

    const VertexMsgMap& GetMsgs() const
    {
        return m_messages;
    }

protected:
    //    long m_id;
    VertexMsgMap m_messages;
    type m_type;
    Vertex m_vertex;

};

class SumNode : public SPGMnode
{
    std::vector<Real> m_weights;
    std::vector<Real> m_logBeta;

public:
    SumNode(Vertex v);
    void SetWeights(const std::vector<Real>& w);
    void SetWeightsRandom();

    virtual void Eval();
    virtual void PassDerivative();

    const std::vector<Real>& GetWeights() const
    {
        return m_weights;
    }

    void StackBeta(const std::vector<Real>& sampleLogProb);
    void PassDerivative_andStackBeta(const std::vector<Real>& sampleLogProb);

    void ResetBeta()
    {
        m_logBeta = std::vector<Real>(m_weights.size());
        for (int i = 0; i < m_logBeta.size(); i++)
            m_logBeta[i] = ZERO_LOGVAL;
    }

    Real GetGradient(int child) const;

    std::string ToStringWeights();

    virtual std::string ToString()
    {
        std::stringstream str;
        str << "+[" << GetVertex() << "]";
        return str.str();
    }

    void EMupdate(Real uniformPrior);

    SumNode* CreateCopy() const
    {
        SumNode* copy = new SumNode(m_vertex);
        copy->m_weights = m_weights;
        CopyData(copy);
        return copy;
    }
};

class ProdNode : public SPGMnode
{
public:
    ProdNode(Vertex v);
    virtual void Eval();
    virtual void PassDerivative();

    virtual std::string ToString()
    {
        std::stringstream str;
        str << "x[" << GetVertex() << "]";
        return str.str();
    }

    ProdNode* CreateCopy() const
    {
        ProdNode* copy = new ProdNode(m_vertex);
        CopyData(copy);
        return copy;
    }
};

class Vnode : public SPGMnode
{
    const int m_var; //index to the var

public:

    Vnode(int var, Vertex v) : SPGMnode(SPGMnode::vnode, v), m_var(var)
    {
        assert(var >= 0);
    };

    int GetVar() const
    {
        return m_var;
    }

    void Eval(const DataMatrix& data);

    Vnode* CreateCopy(const std::map<const CPTfactor*, CPTfactor*>& factorMap) const
    {
        Vnode* copy = new Vnode(m_var, m_vertex);
        CopyData(copy, factorMap);
        return copy;
    }

    void StackBeta(const DataMatrix& data, const std::vector<Real>& sampleLogProb)
    {
        const std::vector<Message*>& childMsgs = m_messages[0]->m_childMsgs; //same children for all output messages
        const std::vector<byte>& dataRow = data.ByCol().GetRow(m_var);
        assert(childMsgs.size() <= 1);
        for (VertexMsgMap::iterator it = m_messages.begin(); it != m_messages.end(); it++)
        {
            Message* msg = *it;
            RealMatrix& logDSpnDmessage = msg->m_logDSpnDmessage;
            CPTfactor* factor = msg->m_factor;

            if (childMsgs.size() == 1)
            {
                factor->StackBeta((childMsgs[0]->m_logMessage), logDSpnDmessage, dataRow, sampleLogProb);
            }
            else //no children
            {
                RealMatrix emptyMat;
                factor->StackBeta(emptyMat, logDSpnDmessage, dataRow, sampleLogProb);
            }

        }
    }

    virtual void PassDerivative()
    {
        assert(false && "todo: remove virtual");
    }
    //    void PassDerivative(const DataMatrix& data);
    void PassDerivative_andStackBeta(const DataMatrix& data, const std::vector<Real>& sampleLogProb);

    void SetFactor(CPTfactor* f, Vertex toVertex)
    {
        if (toVertex == ROOT_VERTEX)
            assert(f->NOut() == 1);
        Message* m = GetMessage(toVertex);
        assert(m != NULL);
        //        assert(f->NOut() == m->m_logMessage.nCols());
        //        assert(f->NOut() == m->m_logDSpnDmessage.nCols());
        m->m_factor = f;
    }

    virtual std::string ToString()
    {
        std::stringstream str;
        str << "Vnode(" << m_var << ")[" << GetVertex() << "]";
        return str.str();
    }
};

