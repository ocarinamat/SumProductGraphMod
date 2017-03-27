//#include "SumMessage.h"
//#include "ProdMessage.h"
//#include "Variable.h"
//
//void SumMessage::Generate(ProdMessage* pMsg)
//{
//    //the output variables of the input message must be a subset of OutV, EXCEPT sumdOut which is not in OutV. Check this.
//    
//    //generate one product node for each assignment to the OutV variables
//    int Nperm = m_perm.NumberOfPermutations();
//
//    //setup one node per each assignment of the output variables
//    m_nodes = new SumNode* [Nperm];
//    for(int i =0; i<Nperm; i++)
//    {
//        m_perm.Ind2perm(i);
//        std::vector<int> outPerm = m_perm.GetSupersetVals();
//        Node node = m_nodes[i];
//        
//        //attach to that node children for every value of the summed out var
//        for (int s=0; i<m_SumdV.m_Nstates; i++) // for every value of the summed out var...
//        {
//            m_perm.SetSVal(s);
//            std::vector<int> inPerm = m_perm.GetSubsetValsPlusS();
//            ProdNode pnode = pMsg->GetNode(inPerm);
//            node.AddChild(pnode);
//        }     
//    }
//}