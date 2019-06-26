/* This software is distributed under the GNU Lesser General Public License */
//==========================================================================
//
//   graph_test.cpp
//
//==========================================================================
// $Id: graph_test.cpp,v 1.1 2003/01/14 16:50:46 raitner Exp $ 

#include <iostream>

#include <GTL/graph.h>
#include <GTL/node_map.h>

#ifdef __GTL_MSVCC
#   ifdef _DEBUG
#	define new DEBUG_NEW
#	undef THIS_FILE
	static char THIS_FILE[] = __FILE__;
#   endif   // _DEBUG
#endif	// __GTL_MSVCC

class vertex {

public:
  vertex() {};
  vertex(int myId) {id=myId;}
  virtual ~vertex() {cout<<" Vertex Detructor ..."<<endl;}
  int id;
};

class parton {

public:
  
  parton() {};
  parton (double mpt, double meta, double mphi, double me, bool mfinal) {pt=mpt;eta=meta;phi=mphi;e=me;final=mfinal;}
  virtual ~parton() {cout<<" Parton Destructor ..."<<endl;}

  bool isFinal() {return final;}
       
  double pt, eta, phi, e;
  bool final;
};

class shower2 : public graph {

public:
  
  shower2() : graph() {cout<<"Shower2 Constrcutor ..."<<endl;}
  virtual ~shower2() {cout<<"Shower2 Detrcutor ..."<<endl;}
  
  node new_vertex(shared_ptr<vertex> v) {node n=graph::new_node();XX[n]=v; return n;}
  //void new_parton(node s, node t, shared_ptr<parton> p) {edge e=graph::new_edge(s,t);PP[e]=p;}//std::move(p);} //; return e;}
  void new_parton(node s, node t, unique_ptr<parton> p) {edge e=graph::new_edge(s,t);PP[e]=std::move(p);}

  //unique_ptr<parton> GetParton(edge e) {return unique_ptr<parton> (PP[e]);}
  
  //int GetNodeValue(node n) const {return XX[n];}
  void save_node_info_handler (ostream *o, node n) const { *o<<"MyID = "<<XX[n]->id<<endl;}
  void save_edge_info_handler (ostream *o, edge n) const { *o<<"Value = "<<PP[n]->pt<<endl;}
  double GetEdgeValue(edge n) const {return PP[n]->pt;}
  int GetNodeValue(node n) const {return XX[n]->id;}
  bool GetFinal(edge e) const {return PP[e]->final;}
  //void pre_new_edge_handler(node s,node t) {};
  //void post_new_edge_handler(edge e) {};
  void pre_clear_handler() {
    edge_iterator git3, gend3;    
    for (git3 = edges_begin(), gend3 = edges_end(); git3 != gend3; ++git3)      
      {
	//cout<<*git3<<" "<<mS->GetEdgeValue(*git3)<<endl;
	PP[*git3]=nullptr;	
      }
    node_iterator git4, gend4; 
    for (git4 = nodes_begin(), gend4 = nodes_end(); git4 != gend4; ++git4)      
      {
	XX[*git4]=nullptr;
      }
  }

  // vector<shared_ptr<vertex>> GetVertices() {};
  /*
  vector<shared_ptr<parton>> GetPartons() {
    vector<shared_ptr<parton>> myv; list<edge> le=all_edges();
    //int n=le.size();
    //for (int i=0;i<n;i++)
    int n=0;

    for (list<edge>::iterator it=le.begin(); it!=le.end(); ++it)
      {
	if (PP[*it]->isFinal())
	  myv.push_back(PP[*it]);
	//cout<<PP[*it]->pt<<endl;
      }
    return myv;
  }
  */
  
  vector<weak_ptr<parton>> GetPartons() {
    vector<weak_ptr<parton>> myv; list<edge> le=all_edges();
    //int n=le.size();
    //for (int i=0;i<n;i++)
    int n=0;

    for (list<edge>::iterator it=le.begin(); it!=le.end(); ++it)
      {
	if (PP[*it]->isFinal())
	  myv.push_back(PP[*it]);
	//cout<<PP[*it]->pt<<endl;
      }
    return myv;
  }
  
 private:
 
  node_map<shared_ptr<vertex>> XX;
  edge_map<shared_ptr<parton>> PP; //unique_ptr not working !???
  
};

int main (int argc, char* argv[])
{
  
    cout << "Loading graph and preserving ids" << endl;
    shower2 G;
    if (G.load("test.gml", true).err_num != GML_OK) {
	cout << "Loading failed" << endl;
	exit(1);
    }

    cout << "Loading OK" << endl;

    if (G.number_of_ids(node()) != 20) {
	cout << "Wrong number of ids: " << G.number_of_ids(node()) << endl;
	exit(1);
    }

    cout << "Number of ids OK" << endl;

    cout << "Loading graph and preserving ids" << endl;
    
    graph G1;
    if (G.load("test.gml", false).err_num != GML_OK) {
	cout << "Loading failed" << endl;
	exit(1);
    }

    cout << "Loading OK" << endl;

    if (G.number_of_ids(node()) != 2) {
	cout << "Wrong number of ids: " << G.number_of_ids(node()) << endl;
	exit(1);
    }

    cout << "Number of ids OK" << endl;
}
