/* This software is distributed under the GNU Lesser General Public License */
//==========================================================================
//
//   bellman_ford_test.cpp
//
//==========================================================================
// $Id: bellman_ford_test.cpp,v 1.2 2002/11/07 13:38:37 raitner Exp $

#include <iostream>

#include <GTL/graph.h>
#include <GTL/bellman_ford.h>
#include <GTL/edge_map.h>
#include <GTL/node_map.h>

#ifdef __GTL_MSVCC
#   ifdef _DEBUG
#	define new DEBUG_NEW
#	undef THIS_FILE
	static char THIS_FILE[] = __FILE__;
#   endif   // _DEBUG
#endif	// __GTL_MSVCC

class my_graph : public graph {

public:

  my_graph() : graph () {}

  node new_vertex(int x) {node n=graph::new_node();X[n]=x; return n;}
  void new_parton(node s, node t, int p) {edge e=graph::new_edge(s,t) ;P[e]=p;} //; return e;}

  void save_node_info_handler (ostream *o, node n) const { *o<<"Value = "<<X[n]<<endl;}
  void save_edge_info_handler (ostream *o, edge n) const { *o<<"Value = "<<P[n]<<endl;}

  int GetNodeValue(node n) const {return X[n];}
  int GetEdgeValue(edge n) const {return P[n];}
  
private:
 
  node_map<int> X;
  edge_map<int> P;
  
};


class parton {

public:
  
  parton() {};
  parton (double mpt, double meta, double mphi, double me, bool mfinal) {pt=mpt;eta=meta;phi=mphi;e=me;final=mfinal;}

  bool isFinal() {return final;}
       
  double pt, eta, phi, e;
  bool final;
};

class shower : public graph {

public:
  
  shower() : graph() {}

  node new_vertex(int x) {node n=graph::new_node();XX[n]=x; return n;}
  void new_parton(node s, node t, parton p) {edge e=graph::new_edge(s,t) ;PP[e]=p;} //; return e;}
  int GetNodeValue(node n) const {return XX[n];}
  void save_edge_info_handler (ostream *o, edge n) const { *o<<"Value = "<<PP[n].pt<<endl;}
  double GetEdgeValue(edge n) const {return PP[n].pt;}
  
 private:
 
  node_map<int> XX;
  edge_map<parton> PP;
  
}; 

int main (int argc, char* argv[])
{
    graph G;
    G.make_directed();
    node n1 = G.new_node();
    node n2 = G.new_node();

    edge e1 = G.new_edge(n1,n2);
    
    
    node_map<int> n(G,1);
    edge_map<int> e(G,10);

    //edge_map<int> e1(n1,n2,10);

    //G.save();

    my_graph S;
    
    node n11=S.new_vertex(0);
    node n22=S.new_vertex(1);
    node n33=S.new_vertex(2);
    node n44=S.new_vertex(3);
    node n55=S.new_vertex(5);

    //edge e11=S.new_parton(n11,n22,10);
    S.new_parton(n11,n22,10);
    S.new_parton(n11,n33,11);
    S.new_parton(n33,n44,20);
    S.new_parton(n33,n55,21);
    
    //S.save();
    
    graph::node_iterator it, end;
    
    for (it = S.nodes_begin(), end = S.nodes_end(); it != end; ++it)      
      {
	cout<<*it<<" "<<S.GetNodeValue((node) *it)<<endl;
      }

    graph::edge_iterator it2, end2;
    
    for (it2 = S.edges_begin(), end2 = S.edges_end(); it2 != end2; ++it2)      
      {
	cout<<*it2<<" "<<S.GetEdgeValue((edge) *it2)<<endl;
      }

    cout<<endl;
    
    shower gS;
    
    node nn11=gS.new_vertex(0);
    node nn22=gS.new_vertex(1);
    node nn33=gS.new_vertex(1);
    node nn44=gS.new_vertex(2);
    node nn55=gS.new_vertex(2);

    parton p(100,0,0,100,false);
    
    gS.new_parton(nn11,nn22,p);
    gS.new_parton(nn11,nn33,parton(200,0,0,200,false));
    gS.new_parton(nn33,nn44,parton(50,0,0,50,true));
    gS.new_parton(nn33,nn55,parton(150,0,0,150,true));
    
    shower::node_iterator git, gend;
    
    for (git = gS.nodes_begin(), gend = gS.nodes_end(); git != gend; ++git)      
      {
	cout<<*git<<" "<<gS.GetNodeValue(*git)<<endl;
      }

    shower::edge_iterator git2, gend2;
    
    for (git2 = gS.edges_begin(), gend2 = gS.edges_end(); git2 != gend2; ++git2)      
      {
	cout<<*git2<<" "<<gS.GetEdgeValue(*git2)<<endl;
      }

    /*
    node n1 = G.new_node();
    node n2 = G.new_node();
    node n3 = G.new_node();
    node n4 = G.new_node();
    node n5 = G.new_node();
    node n6 = G.new_node();

    edge e1_2 = G.new_edge(n1, n2);
    edge e2_3 = G.new_edge(n2, n3);
    edge e1_6 = G.new_edge(n1, n6);
    edge e3_6 = G.new_edge(n3, n6);
    edge e4_3 = G.new_edge(n4, n3);
    edge e6_4 = G.new_edge(n6, n4);
    edge e6_5 = G.new_edge(n6, n5);
    edge e5_4 = G.new_edge(n5, n4);
    edge e5_1 = G.new_edge(n5, n1);

    edge_map<double> w(G);
    w[e1_2] = 1;
    w[e2_3] = 2;
    w[e1_6] = 8;
    w[e3_6] = 3;
    w[e4_3] = 2;
    w[e6_4] = 1;
    w[e6_5] = 3;
    w[e5_4] = 10;
    w[e5_1] = 2;

    node_map<double> result(G);
    result[n1] = 0;
    result[n2] = 1;
    result[n3] = 3;
    result[n4] = 7;
    result[n5] = 9;
    result[n6] = 6;

    node_map<node> preds(G);
    preds[n1] = node();
    preds[n2] = n1;
    preds[n3] = n2;
    preds[n4] = n6;
    preds[n5] = n6;
    preds[n6] = n3;

    bellman_ford B;
    B.store_preds(true);
    B.source(n1);
    B.weights(w);

    G.save("test.gml");
    
    cout << "Bellman Ford with positive edge weights" << endl;
    
    if (B.check(G) == algorithm::GTL_OK) 
    {
	cout << "check OK" << endl; 
    } 
    else 
    {
	cout << "check FAILED" << endl;
	exit(1);
    }

    if (B.run(G) == algorithm::GTL_OK) 
    {
	cout << "run OK" << endl; 
    } 
    else 
    {
	cout << "run FAILED" << endl;
	exit(1);
    }    

    graph::node_iterator it, end;
    
    for (it = G.nodes_begin(), end = G.nodes_end(); it != end; ++it)
    {
	if(result[*it] != B.distance(*it))
	{
	    cout << "Distance for node " << *it << " FAILED" << endl;
	    exit(1);
	}
    }

    cout << "Distances OK" << endl;

    for (it = G.nodes_begin(), end = G.nodes_end(); it != end; ++it)
    {
	if(preds[*it] != B.predecessor_node(*it))
	{
	    cout << "Predecessor for node " << *it << " FAILED" << endl;
	    exit(1);
	}
    }

    cout << "Predecessors OK" << endl;

    if (B.negative_cycle()) 
    {
	cout << "Negative Cycle FAILED" << endl;
    }
    else 
    {
	cout << "Negative Cycle OK" << endl;	
    }
    */
}
