/* 
 * File:   ProdMessage.h
 * Author: mdesana
 *
 * Created on May 5, 2014, 3:04 PM
 */

#ifndef PRODMESSAGE_H
#define	PRODMESSAGE_H

#include "Message.h"
class SumMessage;

//todo delete
class SumNode;

class ProdMessage : Message
{
    VarSet m_indicatorV;
public:
    void Generate(SumMessage* sMsg);
};

#endif	/* PRODMESSAGE_H */

