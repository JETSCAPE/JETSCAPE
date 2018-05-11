#include <iostream>

#include <GTL/graph.h>
#include <GTL/bellman_ford.h>
#include <GTL/edge_map.h>
#include <GTL/node_map.h>
#include <GTL/dfs.h>
#include <memory>
#include <list>
#include <vector>
#include "helper.h"

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
  void save_node_info_handler (ostream *o, node n) const { *o<<"MyID "<<XX[n]->id<<endl;}
  void save_edge_info_handler (ostream *o, edge n) const { *o<<"pT = "<<PP[n]->pt<<endl;}

  void load_edge_info_handler (edge e, GML_pair *read) {
    cout<<"Load edge ... "<<e<<endl;
    struct GML_pair* tmp = read;
    
    while (tmp) {

      //printf ("*KEY* : %s \n", tmp->key);
      if (((string) (tmp->key)).find("pT")<1)
	  {
	  
	    PP[e]=make_shared<parton>(tmp->value.floating,0,0,tmp->value.floating,false);
	  }
	
	tmp = tmp->next;
    }
  }
  
  //maybe smarter to store node id and edge source target in vertex/parton vector and recreate graph from scratch ....
  void load_node_info_handler (node n, GML_pair *read) {
    cout<<"Load node ... "<<n<<endl;
    //cout<<read->key<<endl;//" "<<read->value.floating<<endl;
    struct GML_pair* tmp = read;
    int i;
    int level=2;
    int nC=0;
    
    while (tmp) {
     
      //for (i = 0; i < level; i++) {
      //    printf ("    ");
      //}

	//printf ("*KEY* : %s \n", tmp->key);

	//printf ("%s", tmp->key.id);
	
	//if (((string) (tmp->key)).find("text")) cout<<"text .."<<endl;

	/*
	switch (tmp->kind) {
	case GML_INT:
	    printf ("  *VALUE* (long) : %ld \n", tmp->value.integer);
	    break;

	case GML_DOUBLE:
	    printf ("  *VALUE* (double) : %f \n", tmp->value.floating);
	    break;

	case GML_STRING:
	    printf ("  *VALUE* (string) : %s \n", tmp->value.str);
	    break;
	    
	case GML_LIST:
	    printf ("  *VALUE* (list) : \n");
	    GML_print_list (tmp->value.list, level+1);
	    break;

	default:
	    break;
	}
	*/

	//cout<<((string) (tmp->key)).find("text")<<endl;
	if (((string) (tmp->key)).find("MyID")<1)
	  {
	    //cout<<n<<" "<<nC<<" "<<i<<endl;
	    //cout<<n<<" "<<tmp->value.integer<<endl;
	    XX[n]=make_shared<vertex>(tmp->value.integer);
	  }
	
	//nC++;
	//cout<<nC<<endl;
	tmp = tmp->next;
    }
  }
  
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
    auto mS=make_shared<shower2>();
    //auto G=make_shared<graph>();

    //cout<<G->load("test.gml", true).err_num<<endl;
    if (mS->load("test.gml", true).err_num != GML_OK) {
      cout << "Loading failed" << endl;
	
	//exit(1);
    }

    //G->save();
    //cout << "Loading OK" << endl;

    dfs search;
    search.start_node();// defaulted to first node ... //mS->all_nodes().begin());
    search.run(*mS);
    cout<< search.number_of_reached_nodes () <<endl;

    dfs::dfs_iterator itt2, endt2;
    for (itt2 = search.begin(), endt2=search.end(); itt2 !=endt2; ++itt2)
      {
	cout<<*itt2<<" "<<mS->GetNodeValue(*itt2)<<endl;
      }

    mS=0;
    //mS->clear();
    cout<<"Finished ..."<<endl;
}

