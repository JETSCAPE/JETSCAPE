// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

//Test PartonShower with graphh from GTL
#include "PartonShowerPy8.h"
#include <GTL/dfs.h>
#include <GTL/dijkstra.h>
#include <GTL/node.h>

#include <sstream>
#include <fstream>
#include <iomanip>
#include <random>

#include "MakeUniqueHelper.h"

using namespace std;

namespace Jetscape {

  // quick and dirty fix ...
default_random_engine generator;
uniform_real_distribution<double> distribution(0.2,1.);

PartonShowerPy8::~PartonShowerPy8()
{
  VERBOSESHOWER(8);//pFinal.clear();//pVec.clear();vVec.clear();
}

void PartonShowerPy8::AddTime()
{
  JSDEBUG<<"Add time: Dummy so far; counts generations ...";
  JSINFO<<"  Add time: Dummy so far; counts generations ...";

  if (!isTime) {

    node_iterator nIt,nEnd;
    //vector<node> oneTotwo;
    //vector<node> oneToone;

    //int t=0;

    for (nIt = nodes_begin(), nEnd = nodes_end(); nIt != nEnd; ++nIt)
      {
	node nT=*nIt;

	//cout<<nT<<" "<<GetVertex(nT)->x_in().t()<<" ";
	GTL::node::adj_nodes_iterator aIt,aEnd;

	//double rNum=distribution(generator);

	for (aIt = nT.adj_nodes_begin(), aEnd = nT.adj_nodes_end(); aIt != aEnd; ++aIt)
	  {
	    node nA=*aIt;
	    //cout<<nA<<" ";
	    GetVertex(nA)->x_in().Set(0,0,0,GetVertex(nT)->x_in().t()+1);//+rNum);
	  }
	//cout<<endl;
      }

    //cout<<endl;
    isTime=true;
  }
}

void PartonShowerPy8::CleanPythiaGraph()
{

  JSINFO<<"  Clean Pythia8 Graph of duplicates after transform ...";

  auto eF=GetFinalEdges();
  int nF=eF.size();

  for (auto i : eF) cout<<i.source()<<" ";
  cout<<endl;

  dijkstra  search;
  search.source(GetNodeAt(0));// defaulted to first node ...
  //search.target(eF[0].source());
  search.store_preds(true);
  search.run(*this);

  //cout<<search.source()<<" --> "<<search.target()<<endl;

  for (auto e : eF)
    {
      auto Nitt= search.shortest_path_nodes_begin(e.source());//&GetNodeAt[0]);
      auto Nendt= search.shortest_path_nodes_end(e.source());

      //VERBOSE(7)<<*Nitt<<" "<<*Nendt;

      for (Nitt,Nendt; Nitt !=Nendt; ++Nitt)
	{
	  if (IsOneToOne(*Nitt))
	    {
	      node n=*Nitt;

	      edge eIn=*n.in_edges_begin();
	      edge eOut=*n.out_edges_begin();

	      node nS=eIn.source();
	      node nT=eOut.target();

	      auto p1=GetParton(eIn);
	      auto p2=GetParton(eOut);

	      double deltaPt=abs(p1->pt()-p2->pt());

	      //cout<<*Nitt<<" "<<IsOneToOne(*Nitt)<<" "<<deltaPt<<" ";

	      if (deltaPt<1e-6)
		{
		  //e2=edfs[j];
		  //n2=e2.target();

		  GetVertex(n)=nullptr;
		  hide_node(n);

		  new_parton(nS,nT,p2);

		  //eIn.change_target(eOut.source());

		}
	    }
	}
      //cout<<endl;
    }

    //bfs::tree_edges_iterator Eitt, Eendt;
    //for (Eitt = search.tree_edges_begin(), Eendt=search.tree_edges_end(); Eitt !=Eendt; ++Eitt)
    //{
    //  edge eS=*Eitt;
    //	el.push_back(eS);
    //}

  //CleanPythiaGraphOld();

  /*
  if (!isClean) {

  JSDEBUG<<"Clean Pythia8 Graph of duplicates after transform ...";
  JSDEBUG<<"Use DFS to clean up graph ...";

  INFO<<"  Clean Pythia8 Graph of duplicates after transform ...";

  ClearFinalPartonList();

  }
  */
}

void PartonShowerPy8::CleanPythiaGraphOld()
{
  if (!isClean) {

  JSDEBUG<<"Clean Pythia8 Graph of duplicates after transform ...";
  JSDEBUG<<"Use DFS to clean up graph ...";

  JSINFO<<"  Clean Pythia8 Graph of duplicates after transform ...";

  ClearFinalPartonList();

  dfs search;
  search.calc_comp_num(true);
  search.scan_whole_graph(true);
  search.start_node();// defaulted to first node ...
  search.run(*this);

  // -------------
  // REMARK: Probably not the most efficient way .... !!!!!!!!!!!
  // -------------

  vector<edge> edfs;
  dfs::tree_edges_iterator itt, endt;

  for (itt = search.tree_edges_begin(), endt=search.tree_edges_end(); itt !=endt; ++itt)
    {
      edge eS=*itt;
      edfs.push_back(eS);
    }

  //cout<<edfs.size()-1<<endl;

  for (int i=0;i<edfs.size()-1;i++)
    {
      shared_ptr<Parton> p1=GetParton(edfs[i]);

      if (p1==nullptr) continue;

      int nD=0;

      edge e1=edfs[i];
      node n1=e1.source();

      edge e2; node n2;

      for (int j=i+1;j<edfs.size();j++)
	{

	  //cout<<i<<" "<<j<<endl;
	  //cout<<p1<<" "<<edfs[j]<<endl; //undefined edge ... check for it ...!!!

	  //if (edfs[j].is_hidden()) continue;

	  //cout<<edfs[j]<<endl;

	  shared_ptr<Parton> p2=GetParton(edfs[j]);

	  if (p2==nullptr || p1==nullptr) continue;

	  //cout<<p1<<" "<<p2<<endl;

	  double deltaPt=abs(p1->pt()-p2->pt());

	  //cout<<deltaPt<<endl;

	  if (deltaPt<1e-6)
	    {
	      e2=edfs[j];
	      n2=e2.target();

	      GetVertex(e2.source())=nullptr;
	      del_node(e2.source());

	      nD++;
	    }
	}

      if (nD!=0)
	{
	  new_parton(n1,n2,p1);
	  i+=nD;
	}
    }

  PrintNodes();
  PrintEdges();

  isClean=true;
  }
}

//REMARK: Probably not the same as softdrop since here not 100% sure once flagged if allowed to change once further upstream
// there is a drop !???? Think about it ...
// Check R dropping too ...
void PartonShowerPy8::TransformAndDropPythiaGraph(double zMin, double R)
{
  if (!isTrans) {
    JSDEBUG<<"Transform and Drop wihth zMin = "<<zMin<<" and R = "<<R<<" Pythia8 recursively from the end --> front ...";
    JSINFO<<"  Transform and Drop wihth zMin = "<<zMin<<" and R = "<<R<<" Pythia8 recursively from the end --> front ...";

  // As a first guess, in principle one should redo after trans ...
  double myR=R+0.05;
  shared_ptr<Parton> pSI=GetPartonAt(0);

  //DEBUG:
  //cout<<*pSI<<endl;

  node_iterator nIt,nEnd;
  vector<node> oneTotwo;
  vector<node> oneToone;
  vector<node> dropped;

  for (nIt = nodes_begin(), nEnd = nodes_end(); nIt != nEnd; ++nIt)
    {
      node nT=*nIt;

      if (nT.indeg()==1 && nT.outdeg()==2)
	{
	  oneTotwo.push_back(nT);
	}
      else if (nT.indeg()==1 && nT.outdeg()==1)
	{
	  oneToone.push_back(nT);
	}
    }

  for (int ll=0; ll < oneTotwo.size(); ll++) {

  for (int i=0;i<oneToone.size();i++)
    {
      node nS=oneToone[oneToone.size()-i-1];
      edge eOut=*(nS.out_edges_begin());
      edge eIn=*(nS.in_edges_begin());

      *GetParton(eIn) = *GetParton(eOut);
    }

  for (int i=0;i<oneTotwo.size();i++)
    {
      node nS=oneTotwo[oneTotwo.size()-i-1];
      edge eOut1=*(nS.out_edges_begin());
      edge eOut2=*(++nS.out_edges_begin());
      edge eIn=*(nS.in_edges_begin());

      fjcore::PseudoJet p;
      bool reset=true;

      double pt1=GetParton(eOut1)->GetPseudoJet().pt();
      double pt2=GetParton(eOut2)->GetPseudoJet().pt();

      double dR1=pSI->GetPseudoJet().delta_R(GetParton(eOut1)->GetPseudoJet());
      double dR2=pSI->GetPseudoJet().delta_R(GetParton(eOut2)->GetPseudoJet());

      double z=0;

      if (pt1>pt2)
	z=pt2/(pt1+pt2);
      else
	z=pt1/(pt1+pt2);

      if (dR1<myR && dR2<myR)
	{
	  if (z>zMin)
	    {
	      p = GetParton(eOut1)->GetPseudoJet() + GetParton(eOut2)->GetPseudoJet();
	    }
	  else
	    {
	      //DEBUG:
	      //cout<<z<<" "<<endl;

	      if (pt1>pt2)
		{
		  p=GetParton(eOut1)->GetPseudoJet();
		  if (GetVertex(eOut2.target())->x_in().t()>-99)
		    {
		      GetVertex(eOut2.target())->x_in().Set(-99,-99,-99,-99);
		      dropped.push_back(eOut2.target());
		    }
		}
	      else
		{
		  p=GetParton(eOut2)->GetPseudoJet();
		  if (GetVertex(eOut1.target())->x_in().t()>-99)
		    {
		      GetVertex(eOut1.target())->x_in().Set(-99,-99,-99,-99);
		      dropped.push_back(eOut1.target());
		    }
		}
	    }
	}
      else if (dR1<myR)
	{
	  p=GetParton(eOut1)->GetPseudoJet();
	  if (GetVertex(eOut2.target())->x_in().t()>-99)
	    {
	      GetVertex(eOut2.target())->x_in().Set(-99,-99,-99,-99);
	      dropped.push_back(eOut2.target());
	    }
	}
      else if (dR2<myR)
	{
	  p=GetParton(eOut2)->GetPseudoJet();
	  if (GetVertex(eOut1.target())->x_in().t()>-99)
	    {
	      GetVertex(eOut1.target())->x_in().Set(-99,-99,-99,-99);
	      dropped.push_back(eOut1.target());
	    }
	}
      else
	{
	  reset=false;
	  if (GetVertex(eOut1.target())->x_in().t()>-99)
	    {
	      GetVertex(eOut1.target())->x_in().Set(-99,-99,-99,-99);
	      dropped.push_back(eOut1.target());
	    }

	  if (GetVertex(eOut2.target())->x_in().t()>-99)
	    {
	      GetVertex(eOut2.target())->x_in().Set(-99,-99,-99,-99);
	      dropped.push_back(eOut2.target());
	    }
	}

      if (reset)
	{
	  auto pIn = *GetParton(eIn);

	  // Just 4mom sum ....
	  auto pTrans = make_shared<Parton>(GetParton(eIn)->plabel(),GetParton(eIn)->pid(),0,p.pt(),p.eta(),p.phi(),p.e());

	  *GetParton(eIn) = *pTrans;
	}
    }
  }

  /*
   cout<<"Dropped ..."<<endl;
   for (nIt = nodes_begin(), nEnd = nodes_end(); nIt != nEnd; ++nIt)
     {
       node nT=*nIt;
       cout<<nT<<endl;

       if (GetVertex(nT)->x_in().t()<-98)
	 {
	   //cout<<nT<<endl;
	   //GetVertex(nT)=nullptr;
	   //del_node(nT);
	   //cout<<"del_node ..."<<endl;
	   dropped.push_back(nT);
	 }
     }
   cout<<dropped.size()<<endl;
   cout<<" done"<<endl;
  */

   for (int i=0;i<dropped.size();i++)
     {
       dfs search;
       search.start_node(dropped[i]);
       search.run(*this);

       dfs::dfs_iterator itt2, endt2;
       for (itt2 = search.begin(), endt2=search.end(); itt2 !=endt2; ++itt2)
	 {
	   //DEBUG:
	   //cout<<*itt2<<" ";//<<"="<<search.dfs_num(*itt2)<<" ";
	   hide_node(*itt2);
	 }
     }

   JSINFO<<" Number of partons dropped = "<<dropped.size();

   PrintNodes();
   PrintEdges();

   isTrans=true;
  }
}

// void PartonShowerPy8::TransformAndDropPythiaGraph(double zMin)
// {
//   if (!isTrans) {
//     DEBUG<<"Transform and Drop wihth zMin = "<<zMin<<" Pythia8 recursively from the end --> front ...";

//   node_iterator nIt,nEnd;
//   vector<node> oneTotwo;
//   vector<node> oneToone;
//   vector<node> dropped;

//   for (nIt = nodes_begin(), nEnd = nodes_end(); nIt != nEnd; ++nIt)
//     {
//       node nT=*nIt;

//       if (nT.indeg()==1 && nT.outdeg()==2)
// 	{
// 	  oneTotwo.push_back(nT);
// 	}
//       else if (nT.indeg()==1 && nT.outdeg()==1)
// 	{
// 	  oneToone.push_back(nT);
// 	}
//     }

//   for (int ll=0; ll < oneTotwo.size(); ll++) {

//   for (int i=0;i<oneToone.size();i++)
//     {
//       node nS=oneToone[oneToone.size()-i-1];
//       edge eOut=*(nS.out_edges_begin());
//       edge eIn=*(nS.in_edges_begin());

//       *GetParton(eIn) = *GetParton(eOut);
//     }

//   for (int i=0;i<oneTotwo.size();i++)
//     {
//       node nS=oneTotwo[oneTotwo.size()-i-1];
//       edge eOut1=*(nS.out_edges_begin());
//       edge eOut2=*(++nS.out_edges_begin());
//       edge eIn=*(nS.in_edges_begin());

//       fjcore::PseudoJet p;

//       double pt1=GetParton(eOut1)->GetPseudoJet().pt();
//       double pt2=GetParton(eOut2)->GetPseudoJet().pt();
//       double z=0;

//       if (pt1>pt2)
// 	z=pt2/(pt1+pt2);
//       else
// 	z=pt1/(pt1+pt2);

//       if (z>zMin)
// 	{
// 	  p = GetParton(eOut1)->GetPseudoJet() + GetParton(eOut2)->GetPseudoJet();
// 	}
//       else
// 	{
// 	  //DEBUG:
// 	  //cout<<z<<" "<<endl;

// 	  if (pt1>pt2)
// 	    {
// 	      p=GetParton(eOut1)->GetPseudoJet();
// 	      if (GetVertex(eOut2.target())->x_in().t()>-99)
// 		{
// 		  GetVertex(eOut2.target())->x_in().Set(-99,-99,-99,-99);
// 		  dropped.push_back(eOut2.target());
// 		}
// 	    }
// 	  else
// 	    {
// 	      p=GetParton(eOut2)->GetPseudoJet();
// 	      if (GetVertex(eOut1.target())->x_in().t()>-99)
// 		{
// 		  GetVertex(eOut1.target())->x_in().Set(-99,-99,-99,-99);
// 		  dropped.push_back(eOut1.target());
// 		}
// 	    }
// 	}

//       auto pIn = *GetParton(eIn);

//       // Just 4mom sum ....
//       auto pTrans = make_shared<Parton>(GetParton(eIn)->plabel(),GetParton(eIn)->pid(),0,p.pt(),p.eta(),p.phi(),p.e());

//       *GetParton(eIn) = *pTrans;
//     }
//   }

//   /*
//    cout<<"Dropped ..."<<endl;
//    for (nIt = nodes_begin(), nEnd = nodes_end(); nIt != nEnd; ++nIt)
//      {
//        node nT=*nIt;
//        cout<<nT<<endl;

//        if (GetVertex(nT)->x_in().t()<-98)
// 	 {
// 	   //cout<<nT<<endl;
// 	   //GetVertex(nT)=nullptr;
// 	   //del_node(nT);
// 	   //cout<<"del_node ..."<<endl;
// 	   dropped.push_back(nT);
// 	 }
//      }
//    cout<<dropped.size()<<endl;
//    cout<<" done"<<endl;
//   */

//    for (int i=0;i<dropped.size();i++)
//      {
//        dfs search;
//        search.start_node(dropped[i]);
//        search.run(*this);
//        dfs::dfs_iterator itt2, endt2;
//        for (itt2 = search.begin(), endt2=search.end(); itt2 !=endt2; ++itt2)
// 	 {
// 	   //DEBUG:
// 	   //cout<<*itt2<<" ";//<<"="<<search.dfs_num(*itt2)<<" ";
// 	   hide_node(*itt2);
// 	 }
//      }

//    INFO<<" Number of partons dropped = "<<dropped.size();

//    PrintNodes();
//    PrintEdges();

//    isTrans=true;
//   }
// }


void PartonShowerPy8::TransformPythiaGraph()
{
  TransformPythiaGraphOld();
}

void PartonShowerPy8::TransformPythiaGraphOld()
{
  if (!isTrans) {
  JSDEBUG<<"Transform Pythia8 recursively from the end --> front ...";
  JSINFO<<"  Transform Pythia8 recursively from the end --> front ...";

  node_iterator nIt,nEnd;
  vector<node> oneTotwo;
  vector<node> oneToone;

  for (nIt = nodes_begin(), nEnd = nodes_end(); nIt != nEnd; ++nIt)
    {
      node nT=*nIt;

      if (nT.indeg()==1 && nT.outdeg()==2)
	{
	  oneTotwo.push_back(nT);
	}
      else if (nT.indeg()==1 && nT.outdeg()==1)
	{
	  oneToone.push_back(nT);
	}
    }

  // quick dirty redo as many as 1->2 (not effcient !!!!)
  // -------------
  // REMARK: Probably not the most efficient way .... !!!!!!!!!!!
  // -------------

  for (int ll=0; ll < oneTotwo.size(); ll++) {

  for (int i=0;i<oneToone.size();i++)
    {
      node nS=oneToone[oneToone.size()-i-1];
      edge eOut=*(nS.out_edges_begin());
      edge eIn=*(nS.in_edges_begin());

      *GetParton(eIn) = *GetParton(eOut);
    }

  for (int i=0;i<oneTotwo.size();i++)
    {
      node nS=oneTotwo[oneTotwo.size()-i-1];
      edge eOut1=*(nS.out_edges_begin());
      edge eOut2=*(++nS.out_edges_begin());
      edge eIn=*(nS.in_edges_begin());

      //DEBUG:
      /*
      cout<<  *GetParton(eOut1) << endl;
      cout<< GetParton(eOut1)->m()<< " "<< GetParton(eOut1)->t() <<endl;
      cout<<  *GetParton(eOut2) << endl;
      cout<< GetParton(eOut2)->m()<< " "<< GetParton(eOut2)->t() <<endl;
      */

      auto p = GetParton(eOut1)->GetPseudoJet() + GetParton(eOut2)->GetPseudoJet();
      //auto p = *GetParton(eOut1) + *GetParton(eOut2);
      auto pIn = *GetParton(eIn);

      double newPt=sqrt(pIn.e()*pIn.e()-p.m()*p.m())/cosh(pIn.eta());

      //DEBUG:
      /*
      cout<<" pIn : "<< pIn.pt() << " " << pIn.eta() << " " << pIn.phi() << " "<< pIn.e() << " "<< pIn.m() <<endl;
      cout<<" p   : "<< p.pt() << " " << p.eta() << " " << p.phi() << " "<< p.e() << " "<< p.m() <<endl;
      cout<< (pIn.e()*pIn.e()-p.m()*p.m()) <<endl;
      //cout<<p.pt()<<" "<<newPt<<endl;
      */

      // Just 4mom sum ....
      auto pTrans = make_shared<Parton>(GetParton(eIn)->plabel(),GetParton(eIn)->pid(),GetParton(eIn)->pstat(),p.pt(),p.eta(),p.phi(),p.e());

      /*
      cout<<" ->"<<endl;
      cout<< *pTrans <<endl;
       //  Reset the momentum due to virtuality2
      double velocity[4];

      for(int j=1;j<=3;j++)
      {
	velocity[j] = pTrans->p(j)/pIn.e();//pTrans->e();
      }

      //pTrans->set_jet_v(velocity);

      //double newPl = std::sqrt( (pTrans->e()*pTrans->e() - pTrans->t() )) ;
      double newPl = std::sqrt( (pIn.e()*pIn.e() - pIn.t() )) ;
      //double velocityMod = std::sqrt(std::pow(pTrans->jet_v_.comp(1),2) + std::pow(pTrans->jet_v_.comp(2),2) + std::pow(pTrans->jet_v_.comp(3),2));
      double velocityMod = sqrt(velocity[1]*velocity[1]+velocity[2]*velocity[2]+velocity[3]*velocity[3]);

      newPl = newPl/velocityMod;

      // double newP[4];
      // newP[0] = e();
      // for(int j=1;j<=3;j++) {
      //   newP[j] = newPl*jet_v_.comp(j);
      // }

      pTrans->reset_momentum( newPl*velocity[1], newPl*velocity[2], newPl*velocity[3], pTrans->e() );

      //auto pTrans = make_shared<Parton>(GetParton(eIn)->plabel(),GetParton(eIn)->pid(),0,newPt,p.eta(),p.phi(),p.e());
      */

      //DEBUG:
      /*
      cout<< *GetParton(eIn) <<endl;
      cout<< *pTrans <<endl;
      //cout<< pTrans->form_time()<<endl;
      //cout<< pTrans->m()<<" "<<sqrt(pTrans->t())<<endl;
      cout<<"---"<<endl;
      */

      /*
      //pTrans->set_t(pTrans->m()*pTrans->m());

      cout<< *pTrans <<endl;
      cout<< pTrans->m()<<" "<<sqrt(pTrans->t())<<endl;
      cout<<"---"<<endl;
      */

      *GetParton(eIn) = *pTrans;

      /*
      fjcore::PseudoJet p = *GetParton(eOut1) + *GetParton(eOut2);

      auto pTrans = make_shared<Parton>(GetParton(eIn)->plabel(),GetParton(eIn)->pid(),0,p.pt(),p.eta(),p.phi(),p.e());
      pTrans->reset_PtYPhiM(.pt(),p.rap(),p.phi(),p.m());

      *GetParton(eIn) = *pTrans;
      */

    }

  }

  PrintNodes();
  PrintEdges();

  isTrans=true;
  }
}

//---------------------------------------------------------------------------------------------
//REMARK: Check also for a more generic approach (also for Cloning) EventGraph::FillEvent() !!!
//---------------------------------------------------------------------------------------------

void PartonShowerPy8::FillShower(Pythia8::Event &event, int ps_id)
{
  VERBOSE(7)<<"Fill Pythia Shower for id = "<<ps_id;
  JSINFO<<" Fill Pythia Shower for id = "<<ps_id;

  node vStart;
  node vEnd;

  Pythia8::Particle p = event[ps_id];

  vector<int> dl;
  dl.push_back(ps_id);

  vector<node> vStartVec;
  vector<node> vStartVecTemp;

  vStart=new_vertex(make_shared<Vertex>(0,0,0,0));
  vEnd=new_vertex(make_shared<Vertex>(0,0,0,p.status()));

  new_parton(vStart,vEnd,make_shared<Parton>(ps_id, p.id(),p.status(),p.pT(),p.y(),p.phi(),p.e()));

  vStartVec.push_back(vEnd);

  VERBOSEPARTON(7,*(GetPartonAt(0)));
  //cout<<*(pShower->GetPartonAt(0))<<endl;

  vector<int> dl2;
  vector<int> dl_temp;
  vector<int> sl;

  do {
    for (int k=0;k<dl.size();k++)
      {
	dl2=event[dl[k]].daughterList();

	for (int i=0;i<dl2.size();i++)
	  {
	    //----------------------------------------------------------------
	    //REMARK: To be checked if this is truly "final state" partons!!!
	    //----------------------------------------------------------------
	    if (abs(event[dl2[i]].status())<63)
	      {
		dl_temp.push_back(dl2[i]);
		p=event[dl2[i]];

		vEnd=new_vertex(make_shared<Vertex>(0,0,0,p.status()));
		new_parton(vStartVec[k],vEnd,make_shared<Parton>(event[dl2[i]].index(), p.id(),p.status(),p.pT(),p.y(),p.phi(),p.e()));

		vStartVecTemp.push_back(vEnd);

		if (event[dl2[i]].motherList().size()>1)
		  JSWARN<<"Pythia GTL Shower Fill: More than one mother (to be implemented): "<<event[dl2[i]].motherList().size();
	      }
	  }
	dl2.clear();
      }

    dl.clear();
    dl=dl_temp;
    dl_temp.clear();

    vStartVec.clear();
    vStartVec=vStartVecTemp;
    vStartVecTemp.clear();

  } while (dl.size()>0);

  dl.clear();dl2.clear();dl_temp.clear();

  PrintNodes();
  PrintEdges();
}

}
