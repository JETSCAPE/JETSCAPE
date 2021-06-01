// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

//Test PartonShower with graphh from GTL

#ifndef JETSCAPEGRAPH_H
#define JETSCAPEGRAPH_H

#include "GTL/graph.h"
#include <GTL/edge_map.h>
#include <GTL/node_map.h>
#include "JetClass.h"
//#include "JetScapeParticles.h"
#include "JetScapeLogger.h"
//#include "helper.h"

// REMARK: Still keeping the orginal PartonShower class for now. Even though
//    in the future it should be derived from the templated based class!!!

namespace Jetscape {

  class Vertex;
  class Parton;
  class JetScapeParticleBase;
  class Hadron;

template<class T>
class JetScapeGraph : public graph
{

public:

  JetScapeGraph<T>();
  virtual ~JetScapeGraph<T>();

  node new_vertex(std::shared_ptr<Vertex> v);
  int new_particle(node s, node t, std::shared_ptr<T> p);

  //virtual std::unique_ptr<JetScapeGraph> Clone(); //not yet working for 2--2 !!! Fix it !!!

  void InsertParticle(edge e, std::shared_ptr<Vertex> v, std::shared_ptr<T> p);
  void InsertParticleAfter(edge e, std::shared_ptr<Vertex> v, std::shared_ptr<T> p);
  void InsertEdge(edge e, edge eIns) {};

  std::shared_ptr<Vertex> GetVertex(node n) {return vMap[n];}
  std::shared_ptr<T> GetParticle(edge e) {return pMap[e];}

  std::shared_ptr<T> GetParticleAt(int n);
  std::shared_ptr<Vertex> GetVertexAt(int n);

  edge GetEdge(std::shared_ptr<T> p) {return eMap.at(p);}
  node GetNode(std::shared_ptr<Vertex> v) {return nMap.at(v);}

  node GetNodeAt(int n);
  edge GetEdgeAt(int n);

  bool IsNtoOne(node n);
  bool IsOneToOne(node n);
  bool IsOneToTwo(node n);
  bool IsEndNode(node n);

  bool IsHighSplitEdge(edge e);
  bool IsHighSplitNode(node n);

  double GetSplitDeltaR(node n);
  double GetSplitZ(node n);
  double GetSplitKt(node n);

  double GetNextNodeTime(node n);
  double GetNodeTime(node n); // {return GetVertex(n)->x_in().t();}

  edge GetHighSplitEdge(node n);
  edge GetLowSplitEdge(node n);

  void ReCalculateSplit(node n) {}; //to be implemented (see Py8ShowerPSG.cc for example)

  void GetBfsSortedListOfOneToTwoNodes(vector<node> &nl);
  void GetBfsSortedListOfNodes(vector<node> &nl);
  void GetBfsSortedListOfNodesAndEdges(vector<node> &nl, vector<edge> &el);
  void GetDfsSortedListOfNodes(vector<node> &nl);

  int GetNumberOfParents(int n);
  int GetNumberOfChilds(int n);

  vector<std::shared_ptr<T>> GetFinalParticles();
  vector<fjcore::PseudoJet> GetFinalParticlesForFastJet();

  void ClearFinalParticleList() {pFinal.clear();}

  int GetNumberOfPartons() const {return number_of_edges();}
  int GetNumberOfVertices() const {return number_of_nodes ();}

  void save_node_info_handler (ostream *o, node n) const;
  void save_edge_info_handler (ostream *o, edge n) const;

  void load_edge_info_handler (edge e, GML_pair *read);
  void load_node_info_handler (node n, GML_pair *read);
  void pre_clear_handler();

  void PrintVertices() {PrintNodes(false);}
  void PrintPartons() {PrintEdges(false);}
  void PrintNodes(bool verbose=true);
  void PrintEdges(bool verbose=true);

  void SaveAsGML(string fName) {save(fName.c_str());}
  void SaveAsGV(string fName);
  void SaveAsGraphML(string fName);

private:

  node_map<std::shared_ptr<Vertex>> vMap;
  edge_map<std::shared_ptr<T>> pMap;

  //REMARK: check if needed clearing in destructor ...
  std::map<std::shared_ptr<T>,edge> eMap;
  std::map<std::shared_ptr<Vertex>,node> nMap;

  vector<std::shared_ptr<T>> pFinal;

};

//typedef JetScapeGraph<JetScapeParticleBase> EventGraph;
 typedef JetScapeGraph<Hadron> HadronGraph;

} // end namespace Jetscape
#endif
