// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...


#include "Py8ShowerPSG.h"
#include "PartonShower.h"
#include "JetEnergyLoss.h"
#include "JetScapeLogger.h"
#include <GTL/bfs.h>

#include<iostream>
#include <sstream>

using namespace std;

namespace Jetscape {

void  Py8ShowerPSG::EvolvePartonInTime(double startT, double endT, double deltaT, edge e, shared_ptr<PartonShower> pS, JetEnergyLoss  &j)
{
  //*********************************************************************************************************************************
  //REMARK: Maybe different just take time diff and divide bei deltaT ... so or so, make sure no missing step or double counting ...
  //*********************************************************************************************************************************

  double fractpart, intpart; fractpart = modf (startT , &intpart); double Ts = intpart + ((int) (fractpart*10))/10.;
  double fractpart2, intpart2; fractpart2 = modf (endT , &intpart2); double Te = intpart2 + ((int) (fractpart2*10))/10.;//+deltaT;

  //DEBUG
  //cout<<Ts<<" "<<Te<<" "<<((Te-Ts+deltaT/2.)/deltaT)<<" "<<(int) ((Te-Ts+deltaT/2.)/deltaT)<<endl;

  double currentT=Ts;
  int nTimeSteps=(int) ((Te-Ts+deltaT/2.)/deltaT);

  //DEBUG
  ostringstream timeEvo;
  timeEvo<<"Time evolve : "<<*pS->GetParton(e)<<" : ";
  cout.precision(1);

  vector<Parton> pIn;
  vector<Parton> pOut;

  pIn.push_back(*pS->GetParton(e));

  node nS=e.source();node nE=e.target();

  //Quick test for only checking the first evolution ...
  //DEBUG:
  //if (e.source().id()==0) {

  for (int i=0;i<nTimeSteps;i++)
  {
    currentT = Ts + i*deltaT;

    // --------------------------------------------------

    j.SentInPartons(deltaT,currentT,0,pIn,pOut);

    if (pOut.size()>1) {JSWARN<<"Can not run in conjecture with 1-->2 processes ..."; exit(-1);}

    // --------------------------------------------------

    //DEBUG:
    //timeEvo<<currentT<<" ";

    if (pOut.size()>0)
    {
      //DEBUG:
      //cout<<"pOut : "<<pOut.front()<<endl;

      // --------------------------------------------------
      // Add to Graph ....

      //***************************************************
      // REMARK: If drag like, maybe not add all new nodes
      //  just the whole time integrated change ... (TBD)
      //***************************************************

      //***************************************************
      // Important: Check Position later !!!
      //***************************************************

      auto v=make_shared<Vertex>(0,0,0,currentT);

      pS->InsertParton(e,v,make_shared<Parton>(pOut.front()));
      e=*nE.in_edges_begin ();

      // --------------------------------------------
      // works with Dummy test in DummyEloss ....
      // check again now ...

      if (pIn.size()>1)
      {
        VERBOSESHOWER(7)<<pIn.size()-1<<" new root node(s) to be added ...";

        for (int l=1;l<pIn.size();l++)
        {
          node vNewRootNode=pS->new_vertex(make_shared<Vertex>(0,0,0,currentT-j.GetDeltaT()));

          // ***************************************************
          // Important: Check Position later !!!
          // ***************************************************

          pS->new_parton(vNewRootNode,nE.in_edges_begin ()->source(),make_shared<Parton>(pIn[l]));
          //pS->new_parton(nE.in_edges_begin ()->source(),vNewRootNode,make_shared<Parton>(pIn[l]));
        }
      }

      // --------------------------------------------------

      pIn.clear();
      pIn.push_back(pOut.front());
    }

    pOut.clear();
  }

  //}

  //DEBUG:
  JSDEBUG<<timeEvo.str();
  cout.precision(2);
}

//  ----------------------------------------------------------------------------------------

void Py8ShowerPSG::DoShower(JetEnergyLoss &j)
{
  double tStart=0;
  double currentTime=0;

  JSINFO<<" --> Use Py8ShowerPSG to traverse the transformed Py8 full shower ...";

  //DEBUG:
  cout.precision(2);

  //  ------------------------------------------------------------------------------
  // Implemented Clone/Copy by hand to have full controll !!!! Seems to work!!!
  auto pSTemp=j.GetInitialPartonShower()->Clone();
  shared_ptr<PartonShower> pS=move(pSTemp);

  vector<node> bfslist;
  vector<node> bfslistOrg;

  pS->GetBfsSortedListOfNodes(bfslist);

  // Maybe store some more infos in the Vertex Class!!!
  vector<double> zOrg;
  for (auto i : bfslist) zOrg.push_back(pS->GetSplitZ(i));

  //REMARK: Not the nicest way ... maybe store in Vertex class !???
  int posOrg=0;

  for (auto currNode : bfslist)
  {
    shared_ptr<Vertex> vS=pS->GetVertex(currNode);
    double currTime=vS->x_in().t();

    //cout<<currNode<<" "<<currTime<<endl;

    if (pS->IsOneToTwo(currNode))
    {

      // ------------------------------------------------------------------------------
      // Recalculate split ...
      // Not really correct !!!!! But move on for now ...
      // Move to PartonShower in addition  ... !?

      edge eIn=*currNode.in_edges_begin();
      edge eHigh=pS->GetHighSplitEdge(currNode);
      edge eLow=pS->GetLowSplitEdge(currNode);

      auto pIn=pS->GetParton(eIn);
      auto pH=pS->GetParton(eHigh);
      auto pL=pS->GetParton(eLow);

      double zNew=pS->GetSplitZ(currNode);
      // Not the nicest fix to get original z ...
      double z=zOrg[posOrg];//pSorg->GetSplitZ(bfslistOrg[posOrg]);

      //DEBUG:
      //cout<<currNode<<" "<<posOrg<<" "<<z<<" "<<zNew<<endl;

      // **************************************************
      // REMARK: Not sure if this is how one wants to do it,
      // mass/angle the same, but then new z changes ...
      // TBD how/what is the correct way ...
      // Move on for now, can be easily changed later ...
      // **************************************************

      // not working like this on raw pythia no eloss  ...
      //double mH=pH->m();double mL=pL->m();
      //double eH=sqrt(pIn->pt()*(1-z)*pIn->pt()*(1-z)*cosh(pH->eta())*cosh(pH->eta())+mH*mH);
      //auto pHnew=make_shared<Parton>(pH->plabel(),pH->pid(),pH->pstat(),pIn->pt()*(1-z),pH->eta(),pH->phi(),eH);//,pH->x_in()); (not working ...)

      double mL=sqrt(pIn->e()*pIn->e()-pIn->modp2()); //pL->m();
      double eL=sqrt(pIn->pt()*z*pIn->pt()*z*cosh(pL->eta())*cosh(pL->eta())+mL*mL);
      auto pLnew=make_shared<Parton>(pL->plabel(),pL->pid(),pL->pstat(),pIn->pt()*z,pL->eta(),pL->phi(),eL);

      auto vH=pIn->GetPseudoJet()-pLnew->GetPseudoJet();
      auto pHnew=make_shared<Parton>(pH->plabel(),pH->pid(),pH->pstat(),vH.pt(),vH.eta(),vH.phi(),vH.e());

      // Maybe better with reset_momentum ... anyways, works ...
      *pH=*pHnew; *pL=*pLnew;

      EvolvePartonInTime(currTime,pS->GetNextNodeTime(currNode),j.GetDeltaT(),*currNode.out_edges_begin (),pS,j);
      EvolvePartonInTime(currTime,pS->GetNextNodeTime(currNode),j.GetDeltaT(),*++currNode.out_edges_begin (),pS,j);

      // ------------------------------------------------------------------------------
    }
    else if (j.GetInitialPartonShower()->IsEndNode(currNode))
    {
      //****************
      // To be done ...
      //****************
      //DEBUG:
      JSDEBUG<<"1 -> END = "<<currNode<<" evolve until Tmax ...";

    }
    else
    {
      //DEBUG:
      JSDEBUG<<"1 -> 1 = "<<currNode;
      EvolvePartonInTime(currTime,pS->GetNextNodeTime(currNode),j.GetDeltaT(),*currNode.out_edges_begin (),pS,j);
    }

    posOrg++;
  }

  //  ------------------------------------------------------------------------------

  //Attached modifed shower copy to Eloss ...
  j.SetShower(pS);
}

} // end namespace Jetscape
