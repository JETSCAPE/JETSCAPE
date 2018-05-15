/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 * 
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

//PartonShower with graph from GTL

#ifndef PARTONSHOWER_H
#define PARTONSHOWER_H

#include "GTL/graph.h"
#include <GTL/edge_map.h>
#include <GTL/node_map.h>
#include "JetClass.h"
#include "JetScapeLogger.h"

using std::shared_ptr;

namespace Jetscape {

// Think about best interface and what is truly needed, maybe even better
// if no graph at all write a converter function/class and split parton in base
// and after transformer class. TBD ...
  class Vertex;
  class Parton;
  
class PartonShower : public graph
{

public:

  PartonShower();
  virtual ~PartonShower();

  node new_vertex(shared_ptr<Vertex> v);
  int new_parton(node s, node t, shared_ptr<Parton> p);
  
  shared_ptr<Vertex> GetVertex(node n) {return vMap[n];}
  shared_ptr<Parton> GetParton(edge e) {return pMap[e];}

  shared_ptr<Parton> GetPartonAt(int n);
  shared_ptr<Vertex> GetVertexAt(int n);

  node GetNodeAt(int n);
  edge GetEdgeAt(int n);
  
  int GetNumberOfParents(int n);
  int GetNumberOfChilds(int n);

  vector<shared_ptr<Parton>> GetFinalPartons();
  vector<fjcore::PseudoJet> GetFinalPartonsForFastJet();
  
  //unique_ptr<Parton> GetPartonAt(int i);
  //unique_ptr<Vertex> GetVertexAt(int i);
  //vector<unique_ptr<Parton>> GetPartons() {};

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
  
  node_map<shared_ptr<Vertex>> vMap;
  edge_map<shared_ptr<Parton>> pMap;

  vector<shared_ptr<Parton>> pFinal;
  
  //Check map data format (pointer to edge/node !??)
  //map<weak_ptr<Parton>, edge> pToEdgeMap;
  //map<weak_ptr<Vertex>, node> vToNodeMap;

  //void CreateMaps();
  //bool mapsFilled;
  
  //Check here memory issues to avoid cyclic references (use unique pointers instead for return functions!?
  //In general rethink and clean up pointer types for efficiency and safety ...
  //Only fill when needed via Getters ...
  //Can also be done via lists, a bit slower (interface question ...) 
  //vector<weak_ptr<Parton>> pVec;
  //vector<weak_ptr<Vertex>> vVec;

  //void FillVertexVec();
  //void FillPartonVec();
  
};

} // end namespace Jetscape
#endif
