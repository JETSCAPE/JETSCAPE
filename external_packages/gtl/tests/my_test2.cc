#include <iostream>

#include <GTL/graph.h>
#include <GTL/bellman_ford.h>
#include <GTL/edge_map.h>
#include <GTL/node_map.h>
#include <GTL/dfs.h>
#include <memory>
#include <list>

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
  ~parton() {cout<<" Parton Destructor ..."<<endl;}

  bool isFinal() {return final;}
       
  double pt, eta, phi, e;
  bool final;
};

class shower : public graph {

public:
  
  shower() : graph() {cout<<"Shower Constrcutor ..."<<endl;}
  ~shower() {cout<<"Shower Detrcutor ..."<<endl;}
  
  node new_vertex(int x) {node n=graph::new_node();XX[n]=x; return n;}
  void new_parton(node s, node t, parton p) {edge e=graph::new_edge(s,t) ;PP[e]=p;} //; return e;}
  int GetNodeValue(node n) const {return XX[n];}
  void save_edge_info_handler (ostream *o, edge n) const { *o<<"Value = "<<PP[n].pt<<endl;}
  double GetEdgeValue(edge n) const {return PP[n].pt;}
  
 private:
 
  node_map<int> XX;
  edge_map<parton> PP;
  
};

class shower2 : public graph {

public:
  
  shower2() : graph() {cout<<"Shower2 Constrcutor ..."<<endl;}
  virtual ~shower2() {cout<<"Shower2 Detrcutor ..."<<endl;}
  
  node new_vertex(int x) {node n=graph::new_node();XX[n]=x; return n;}
  //void new_parton(node s, node t, shared_ptr<parton> p) {edge e=graph::new_edge(s,t);PP[e]=p;}//std::move(p);} //; return e;}
  void new_parton(node s, node t, unique_ptr<parton> p) {edge e=graph::new_edge(s,t);PP[e]=std::move(p);}
  
  //int GetNodeValue(node n) const {return XX[n];}
  void save_edge_info_handler (ostream *o, edge n) const { *o<<"Value = "<<PP[n]->pt<<endl;}
  double GetEdgeValue(edge n) const {return PP[n]->pt;}
  int GetNodeValue(node n) const {return XX[n];}
  bool GetFinal(edge e) const {return PP[e]->final;}
  //void pre_new_edge_handler(node s,node t) {};
  void post_new_edge_handler(edge e) {};
  void pre_clear_handler() {
    edge_iterator git3, gend3;    
    for (git3 = edges_begin(), gend3 = edges_end(); git3 != gend3; ++git3)      
      {
	//cout<<*git3<<" "<<mS->GetEdgeValue(*git3)<<endl;
	PP[*git3]=nullptr;
      }
  }
  
 private:
 
  node_map<int> XX;
  edge_map<shared_ptr<parton>> PP; //unique_ptr not working !???
  
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
    //S.insert_reverse_edges();
    
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

    cout<<endl;
    
    node nn11=gS.new_vertex(0);
    node nn22=gS.new_vertex(1);
    node nn33=gS.new_vertex(1);
    node nn44=gS.new_vertex(2);
    node nn55=gS.new_vertex(2);

    //parton p(100,0,0,100,false);
    
    gS.new_parton(nn11,nn22,parton(100,0,0,100,false));
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

    //shower *mS=new shower();
    cout<<endl;
    auto mS=make_shared<shower2>();
    mS->make_directed();
    
    node nnn11=mS->new_vertex(0);
    node nnn22=mS->new_vertex(1);
    node nnn33=mS->new_vertex(1);
    node nnn44=mS->new_vertex(2);
    node nnn55=mS->new_vertex(2);
    node nnn66=mS->new_vertex(3);
    node nnn77=mS->new_vertex(3);
    node nnn88=mS->new_vertex(4);
    
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
	cout<<*git3<<" "<<mS->GetEdgeValue(*git3)<<endl;
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
    
    mS->clear(); //does not call parton destructor but with handler it does  ...
    cout<<" mS->clear()"<<endl;
    //mS->save();
    
    mS=0;
    cout<<"Done .."<<endl;
}
