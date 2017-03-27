//#include "SumMessage.h"
//#include "ProdMessage.h"
//
//void ProdMessage::Generate(SumMessage* sMsg)
//{
//    //check that the input messages have disjoint OutV scopes
//    
//    int Nperm = m_perm.NumberOfPermutations();
//
//    //setup one node per each assignment of the output variables
//    m_nodes = new SumNode* [Nperm];
//    for(int i =0; i<Nperm; i++)
//    {
//        m_perm.Ind2perm(i);
//        outPerm = m_perm.GetAVals();
//        node = m_nodes[i];
//        
//        //attach to that node children for every value of the summed out var
//        for (int j=0; j<NinputMessages; j++) // for every value of the summed out var...
//        {
//            inPerm = m_perm.GetSubsetVals(j);
//            ProdNode pnode = sMsg[i]->GetNode(inPerm);
//            node.AddChild(pnode);
//        }     
//    }
//}
