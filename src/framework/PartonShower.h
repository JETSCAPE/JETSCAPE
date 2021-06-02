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
//#include "JetScapeParticles.hpp"
#include "JetScapeLogger.h"
//#include "helper.h"

#include <map>

using std::shared_ptr;

namespace Jetscape {

  class Vertex;
  class Parton;

class PartonShower : public graph
{

public:

  PartonShower();
  virtual ~PartonShower();

  node new_vertex(std::shared_ptr<Vertex> v);
  int new_parton(node s, node t, std::shared_ptr<Parton> p);

  virtual std::unique_ptr<PartonShower> Clone(); //not yet working for 2--2 !!! Fix it !!!

  void InsertParton(edge e, std::shared_ptr<Vertex> v, std::shared_ptr<Parton> p);
  void InsertPartonAfter(edge e, std::shared_ptr<Vertex> v, std::shared_ptr<Parton> p);
  void InsertEdge(edge e, edge eIns) {};

  std::shared_ptr<Vertex> GetVertex(node n) {return vMap[n];}
  std::shared_ptr<Parton> GetParton(edge e) {return pMap[e];}

  std::shared_ptr<Parton> GetPartonAt(int n);
  std::shared_ptr<Vertex> GetVertexAt(int n);

  edge GetEdge(std::shared_ptr<Parton> p) {return eMap.at(p);}
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

  vector<std::shared_ptr<Parton>> GetFinalPartons();
  vector<fjcore::PseudoJet> GetFinalPartonsForFastJet();

  vector<std::shared_ptr<Parton>> GetPartons();

  vector<edge> GetFinalEdges();

  void ClearFinalPartonList() {pFinal.clear();}
  void ClearPartonList() {pAll.clear();}

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
  edge_map<std::shared_ptr<Parton>> pMap;

  std::map<std::shared_ptr<Parton>,edge> eMap;
  std::map<std::shared_ptr<Vertex>,node> nMap;

  vector<std::shared_ptr<Parton>> pFinal;
  vector<std::shared_ptr<Parton>> pAll;

  vector<edge> eFinal;
};

} // end namespace Jetscape
#endif
