#include <iostream>

#include <GTL/graph.h>
#include <GTL/bellman_ford.h>
#include <GTL/edge_map.h>
#include <GTL/node_map.h>
#include <GTL/dfs.h>
#include <GTL/bfs.h>
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
  void save_edge_info_handler (ostream *o, edge n) const { *o<<"pT "<<PP[n]->pt<<endl;}

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
    auto mS=make_shared<shower2>();
    mS->make_directed();

    node nnn11=mS->new_vertex(make_shared<vertex>(0)); //if c++14 or c++11 and helper.h ...
    //node nnn11=mS->new_vertex(unique_ptr<vertex>(new vertex(0)));
    node nnn22=mS->new_vertex(make_shared<vertex>(1));
    node nnn33=mS->new_vertex(make_shared<vertex>(1));
    node nnn44=mS->new_vertex(make_shared<vertex>(2));
    node nnn55=mS->new_vertex(make_shared<vertex>(2));
    node nnn66=mS->new_vertex(make_shared<vertex>(3));
    node nnn77=mS->new_vertex(make_shared<vertex>(3));
    node nnn88=mS->new_vertex(make_shared<vertex>(4));
    
    /*
    node nnn11=mS->new_vertex(make_unique<vertex>(1)); //if c++14 or c++11 and helper.h ...
    //node nnn11=mS->new_vertex(unique_ptr<vertex>(new vertex(0)));
    node nnn22=mS->new_vertex(unique_ptr<vertex>(new vertex(1)));
    node nnn33=mS->new_vertex(unique_ptr<vertex>(new vertex(1)));
    node nnn44=mS->new_vertex(unique_ptr<vertex>(new vertex(2)));
    node nnn55=mS->new_vertex(unique_ptr<vertex>(new vertex(2)));
    node nnn66=mS->new_vertex(unique_ptr<vertex>(new vertex(3)));
    node nnn77=mS->new_vertex(unique_ptr<vertex>(new vertex(3)));
    node nnn88=mS->new_vertex(unique_ptr<vertex>(new vertex(4)));
    */
    
    mS->new_parton(nnn11,nnn22,unique_ptr<parton>(new parton(100,0,0,100,false)));
    mS->new_parton(nnn11,nnn33,unique_ptr<parton>(new parton(200,0,0,200,false)));
    mS->new_parton(nnn33,nnn44,unique_ptr<parton>(new parton(50,0,0,50,true)));
    mS->new_parton(nnn33,nnn55,unique_ptr<parton>(new parton(150,0,0,150,false)));
    mS->new_parton(nnn55,nnn66,unique_ptr<parton>(new parton(80,0,0,80,true)));
    mS->new_parton(nnn55,nnn77,unique_ptr<parton>(new parton(70,0,0,70,false)));
    mS->new_parton(nnn77,nnn88,unique_ptr<parton>(new parton(60,0,0,60,true)));

    //reverse changes all counts --> new inverse graph ... ? (shared pointers ...)
    /*
    mS->new_parton(nnn22,nnn11,unique_ptr<parton>(new parton(100,0,0,100,false)));
    mS->new_parton(nnn33,nnn11,unique_ptr<parton>(new parton(200,0,0,200,false)));
    mS->new_parton(nnn44,nnn33,unique_ptr<parton>(new parton(50,0,0,50,true)));
    mS->new_parton(nnn55,nnn33,unique_ptr<parton>(new parton(150,0,0,150,false)));
    mS->new_parton(nnn66,nnn55,unique_ptr<parton>(new parton(80,0,0,80,true)));
    mS->new_parton(nnn77,nnn55,unique_ptr<parton>(new parton(70,0,0,70,true)));
    */
    
    //mS->new_parton(nnn11,nnn22,make_shared<parton>(100,0,0,100,false));
    //mS->save();

    //mS->insert_reverse_edges(); //problem with new parton creation .... so post/pre handler
    // not really doing the trick
    
    shower2::node_iterator git0, gend0;
    
    for (git0 = mS->nodes_begin(), gend0 = mS->nodes_end(); git0 != gend0; ++git0)      
      {
	cout<<*git0<<" "<<mS->GetNodeValue(*git0)<<" "<<git0->indeg()<<" "<<git0->outdeg()<<endl;
      }
    
    shower2::edge_iterator git3, gend3;
    
    for (git3 = mS->edges_begin(), gend3 = mS->edges_end(); git3 != gend3; ++git3)      
      {
	cout<<*git3<<" "<<mS->GetEdgeValue(*git3)<<" "<<mS->GetFinal(*git3)<<endl;
      }

    list<node> ln=mS->all_nodes();
    cout<<ln.size()<<endl;

    list<edge> le=mS->all_edges(); 
    cout<<le.size()<<endl;
    //le.remove_if( // hmmm, not with parton in graph (inheritence !??) Check ..

    dfs search;
    search.start_node(nnn11);
    search.run(*mS);
    cout<< search.number_of_reached_nodes () <<endl;

    dfs::dfs_iterator itt2, endt2;
    for (itt2 = search.begin(), endt2=search.end(); itt2 !=endt2; ++itt2)
      {
	cout<<*itt2<<endl;
      }
    
    dfs::tree_edges_iterator itt, endt;
    for (itt = search.tree_edges_begin(), endt=search.tree_edges_end(); itt !=endt; ++itt)
      {
	cout<<*itt<<endl;
      }

    bfs search2;
    search2.start_node(nnn11);
    search2.run(*mS);
    cout<< search2.number_of_reached_nodes () <<endl;

    bfs::bfs_iterator itt22, endt22;
    for (itt22 = search2.begin(), endt22=search2.end(); itt22 !=endt22; ++itt22)
      {
	cout<<*itt22<<endl;
      }
    
    // -----
    // why only with shared vectors ? Not reall happy!?? explore with unique_ptr !???
    // somehow issues with the GTL !???

    /*
    vector<shared_ptr<parton>> myv=mS->GetPartons(); // this is dangerous with shared pointers ... !!!!
    cout<<myv.size()<<endl; 
    myv.clear();
    */
    
    vector<weak_ptr<parton>> myv=mS->GetPartons(); //better with weak ... !!!!
    cout<<myv.size()<<endl;
    //if (myv.size()>0)
    for (const auto& v : myv)
      cout<<v.lock()->pt<<endl;
      //cout<<myv[0].lock()->pt<<endl;
    //myv.clear();

    mS->save("test.gml");
    
    mS->clear(); //does not call parton destructor but with handler it does  ...
    cout<<" mS->clear()"<<endl; 
    //mS->save();
    
    mS=0;

    // --------------
    // not really working with inhereting edge and node ...
    /*
    cout<<"  -------------- "<<endl;
    
    auto mS3=make_shared<shower3>();
    mS3->make_directed();

    mS3->new_vertex(0);
    
    mS3->save();
    
    cout<<"Done .."<<endl;
    */
}

/*
// not really working though ... !????
class vertex : public node
{

public :
  
  vertex() : node() {};
  vertex(int mylabel) : node() {label=mylabel;}
  virtual ~vertex() {cout<<"Vertex destructor ..."<<endl;}
  
  int label;
};

class parton3 : public edge
{

public:
  
  parton3() : edge() {};
  parton3 (double mpt, double meta, double mphi, double me, bool mfinal) :edge () {pt=mpt;eta=meta;phi=mphi;e=me;final=mfinal;}
  virtual ~parton3() {cout<<" Parton3 Destructor ..."<<endl;}

  bool isFinal() {return final;}
       
  double pt, eta, phi, e;
  bool final;
};
*/
/*
class shower3 : public graph {

public:
  
  shower3() : graph() {cout<<"Shower3 Constrcutor ..."<<endl;}
  virtual ~shower3() {cout<<"Shower3 Detrcutor ..."<<endl;}
  
  //node new_vertex(int x) {node n=graph::new_node();XX[n]=x; return n;}
  //void new_parton(node s, node t, unique_ptr<parton> p) {edge e=graph::new_edge(s,t);PP[e]=std::move(p);}
  //void save_edge_info_handler (ostream *o, edge n) const { *o<<"Value = "<<PP[n]->pt<<endl;}

  //vertex new_vertex(int i) {unique_ptr<vertex> n = unique_ptr<vertex>(new vertex(i)); *(static_cast<edge>(n)) = graph::new_node();return n;} //XX[n]=x; return n;}

  void new_vertex(int i) {shared_ptr<vertex> n = make_shared<vertex>(i); (node) (*n) = graph::new_node();}//return n;} //XX[n]=x; return n;}
  void new_parton(node s, node t, unique_ptr<parton> p) {edge e=graph::new_edge(s,t);} ;//PP[e]=std::move(p);}
  /*
  //void save_edge_info_handler (ostream *o, edge n) const { *o<<"Value = "<<PP[n]->pt<<endl;}
  
  //double GetEdgeValue(edge n) const {return PP[n]->pt;}
  //int GetNodeValue(node n) const {return XX[n];}
  //bool GetFinal(edge e) const {return PP[e]->final;}
  //void pre_new_edge_handler(node s,node t) {};
  //void post_new_edge_handler(edge e) {};

  /*
  void pre_clear_handler() {
    edge_iterator git3, gend3;    
    for (git3 = edges_begin(), gend3 = edges_end(); git3 != gend3; ++git3)      
      {
	//cout<<*git3<<" "<<mS->GetEdgeValue(*git3)<<endl;
	PP[*git3]=nullptr;
      }
  }
  */
/*  
 private:
 
  //node_map<int> XX;
  //edge_map<shared_ptr<parton>> PP; //unique_ptr not working !???
  
};
*/
