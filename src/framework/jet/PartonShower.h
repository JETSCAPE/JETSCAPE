// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

//Test PartonShower with graphh from GTL

#ifndef PARTONSHOWER_H
#define PARTONSHOWER_H

#include "GTL/graph.h"
#include <GTL/edge_map.h>
#include <GTL/node_map.h>
#include "JetClass.hpp"
#include "JetScapeLogger.h"

class PartonShower : public graph
{

public:

  PartonShower() : graph() {VERBOSESHOWER(8);}
  virtual ~PartonShower() {VERBOSESHOWER(8);}

  node new_vertex(shared_ptr<VertexBase> v);
  void new_parton(node s, node t, shared_ptr<Parton> p);
  weak_ptr<VertexBase> GetVertex(node n) {return vMap[n];}
  weak_ptr<Parton> GetParton(edge e) {return pMap[e];}

  void save_node_info_handler (ostream *o, node n) const; 
  void save_edge_info_handler (ostream *o, edge n) const;

  void load_edge_info_handler (edge e, GML_pair *read);
  void load_node_info_handler (node n, GML_pair *read);
  void pre_clear_handler();

  void PrintVertices() {};
  void PrintPartons() {};
  void PrintNodes();
  void PrintEdges();
  
  void SaveAsGML(string fName) {save(fName.c_str());}
  void SaveAsGV(string fName);
		
private:

  node_map<shared_ptr<VertexBase>> vMap;
  edge_map<shared_ptr<Parton>> pMap; 
  
};
#endif
