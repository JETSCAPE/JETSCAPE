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

#include "PartonShower.h"
#include <sstream>
#include <fstream>
#include <iomanip>
#include "MakeUniqueHelper.h"

using std::setprecision;
using std::fixed;
using std::to_string;
using std::ostringstream;
using std::stringstream;

namespace Jetscape {

PartonShower::PartonShower() : graph() { VERBOSESHOWER(8); }

node PartonShower::new_vertex(shared_ptr<Vertex> v) {
  node n = graph::new_node();
  vMap[n] = v;
  return n;
}

int PartonShower::new_parton(node s, node t, shared_ptr<Parton> p) {
  edge e = graph::new_edge(s, t);
  pMap[e] = p;
  return e.id();
}

/*
void PartonShower::FillPartonVec()
{
  VERBOSE(8);
  
  edge_iterator eIt, eEnd;
  
  for (eIt = edges_begin(), eEnd = edges_end(); eIt != eEnd; ++eIt)
    pVec.push_back(pMap[*eIt]);     
}

void PartonShower::FillVertexVec()
{
  VERBOSE(8);
  
  node_iterator nIt, nEnd;
  
  for (nIt = nodes_begin(), nEnd = nodes_end(); nIt != nEnd; ++nIt)    
    vVec.push_back(vMap[*nIt]);  
}
*/

/*
shared_ptr<Parton> PartonShower::GetPartonAt(int i)
{
  if (pVec.size()==0)
    FillPartonVec();  
  
  return pVec[i].lock();
}

shared_ptr<Vertex> PartonShower::GetVertexAt(int i)
{
  if (vVec.size()==0)
    FillVertexVec();

  return vVec[i].lock();
}
*/

/*
unique_ptr<Parton> PartonShower::GetPartonAt(int i)
{
  if (pVec.size()==0)
    FillPartonVec();  
  
  return make_unique<Parton>(*(pVec[i].lock()));
}

unique_ptr<Vertex> PartonShower::GetVertexAt(int i)
{
  if (vVec.size()==0)
    FillVertexVec();

  return make_unique<Vertex>(*(vVec[i].lock()));
}
*/
/*
void PartonShower::CreateMaps()
{
  VERBOSESHOWER(8);
  edge_iterator eIt, eEnd;
  for (eIt = edges_begin(), eEnd = edges_end(); eIt != eEnd; ++eIt)
    {
      pToEdgeMap[pMap[*eIt]]=*eIt;
    }
  
  node_iterator nIt, nEnd;
  for (nIt = nodes_begin(), nEnd = nodes_end(); nIt != nEnd; ++nIt)
    {
      //vToNodeMap[vMap[*nIt]]=*nIt;
    }
}
*/

vector<shared_ptr<Parton>> PartonShower::GetFinalPartons() {
  //VERBOSESHOWER(8)<<pFinal.size();
  if (pFinal.size() == 0) {
    edge_iterator eIt, eEnd;
    for (eIt = edges_begin(), eEnd = edges_end(); eIt != eEnd; ++eIt) {
      if (eIt->target().outdeg() < 1) {
        if (pMap[*eIt]->pstat() > -10) {
          pFinal.push_back(pMap[*eIt]);
        }
        // DEBUG
        //cout<<eIt->target()<<endl;
      }
    }
    return pFinal;
  } else {
    return pFinal;
  }
}

vector<fjcore::PseudoJet> PartonShower::GetFinalPartonsForFastJet() {
  vector<shared_ptr<Parton>> mP =
      GetFinalPartons(); // to ensure that pFinal is filled

  vector<fjcore::PseudoJet> forFJ;

  for (int i = 0; i < mP.size(); i++)
    forFJ.push_back((*mP[i]).GetPseudoJet());

  return forFJ;
}

int PartonShower::GetNumberOfParents(int n) {
  return GetEdgeAt(n).source().indeg();
}

int PartonShower::GetNumberOfChilds(int n) {
  return GetEdgeAt(n).target().outdeg();
}

edge PartonShower::GetEdgeAt(int n) {
  edge_iterator eIt;
  eIt = edges_begin();
  advance(eIt, n);
  return *eIt;
}

node PartonShower::GetNodeAt(int n) {
  node_iterator nIt;
  nIt = nodes_begin();
  advance(nIt, n);
  return *nIt;
}

shared_ptr<Parton> PartonShower::GetPartonAt(int n) {
  edge_iterator eIt;
  eIt = edges_begin();
  advance(eIt, n);
  return GetParton(*eIt);
}

shared_ptr<Vertex> PartonShower::GetVertexAt(int n) {
  node_iterator nIt;
  nIt = nodes_begin();
  advance(nIt, n);
  return GetVertex(*nIt);
}

PartonShower::~PartonShower() {
  VERBOSESHOWER(8);
  pFinal.clear(); //pVec.clear();vVec.clear();
}

void PartonShower::save_node_info_handler(ostream *o, node n) const {
  *o << "label "
     << "\"" << n.id() << "(" << fixed << setprecision(2) << vMap[n]->x_in().t()
     << ")\"" << endl;
  *o << "x " << vMap[n]->x_in().x() << endl;
  *o << "y " << vMap[n]->x_in().y() << endl;
  *o << "z " << vMap[n]->x_in().z() << endl;
  *o << "t " << vMap[n]->x_in().t() << endl;
}

void PartonShower::save_edge_info_handler(ostream *o, edge e) const {
  *o << "label "
     << "\"(" << fixed << setprecision(2) << pMap[e]->pt() << ")\"" << endl;
  *o << "plabel " << pMap[e]->plabel() << endl;
  *o << "pid " << pMap[e]->pid() << endl;
  *o << "pstat " << pMap[e]->pstat() << endl;
  *o << "pT " << pMap[e]->pt() << endl;
  *o << "eta " << pMap[e]->eta() << endl;
  *o << "phi " << pMap[e]->phi() << endl;
  *o << "E " << pMap[e]->e() << endl;
}

void PartonShower::pre_clear_handler() {
  VERBOSESHOWER(8);
  edge_iterator eIt, eEnd;
  for (eIt = edges_begin(), eEnd = edges_end(); eIt != eEnd; ++eIt) {
    pMap[*eIt] = nullptr;
  }

  node_iterator nIt, nEnd;
  for (nIt = nodes_begin(), nEnd = nodes_end(); nIt != nEnd; ++nIt) {
    vMap[*nIt] = nullptr;
  }
}

void PartonShower::PrintNodes(bool verbose) {
  node_iterator nIt, nEnd;
  ostringstream os;

  if (verbose && JetScapeLogger::Instance()->GetVerboseLevel() > 8) {
    for (nIt = nodes_begin(), nEnd = nodes_end(); nIt != nEnd; ++nIt)
      os << *nIt << "=" << vMap[*nIt]->x_in().t() << " ";
    VERBOSESHOWER(8) << os.str();
  }

  if (!verbose) {
    for (nIt = nodes_begin(), nEnd = nodes_end(); nIt != nEnd; ++nIt)
      os << *nIt << "=" << vMap[*nIt]->x_in().t() << " ";
    cout << "Vertex list : " << os.str() << endl;
  }
}

void PartonShower::PrintEdges(bool verbose) {
  edge_iterator eIt, eEnd;
  ostringstream os;

  if (verbose && JetScapeLogger::Instance()->GetVerboseLevel() > 8) {
    for (eIt = edges_begin(), eEnd = edges_end(); eIt != eEnd; ++eIt)
      os << *eIt << "=" << pMap[*eIt]->pt() << " ";
    VERBOSESHOWER(8) << os.str();
  }

  if (!verbose) {
    for (eIt = edges_begin(), eEnd = edges_end(); eIt != eEnd; ++eIt)
      os << *eIt << "=" << pMap[*eIt]->pt() << " ";
    cout << "Parton list : " << os.str() << endl;
  }
}

// To be extended to store all infos like this in GML
// or store vectore of partons and vertecies with relevant graph infos
// to construct the graph later .... (TBD)
// These handlers are needed if one wants to use load and GML files ...
// Obsolete with new JetScape reader ...

void PartonShower::load_edge_info_handler(edge e, GML_pair *read) {
  VERBOSESHOWER(8) << "Load edge ... " << e;
  struct GML_pair *tmp = read;

  double pT, eta, phi, E;
  int plabel, pid, pstat;

  while (tmp) {
    //printf ("*KEY* : %s \n", tmp->key);
    if (((string)(tmp->key)).find("pT") < 1)
      pT = tmp->value.floating;
    if (((string)(tmp->key)).find("eta") < 1)
      eta = tmp->value.floating;
    if (((string)(tmp->key)).find("phi") < 1)
      phi = tmp->value.floating;
    if (((string)(tmp->key)).find("E") < 1)
      E = tmp->value.floating;
    if (((string)(tmp->key)).find("plabel") < 1)
      plabel = tmp->value.integer;
    if (((string)(tmp->key)).find("pid") < 1)
      pid = tmp->value.integer;
    if (((string)(tmp->key)).find("pstat") < 1)
      pstat = tmp->value.integer;

    tmp = tmp->next;
  }

  pMap[e] = make_shared<Parton>(plabel, pid, pstat, pT, eta, phi, E);
}

void PartonShower::load_node_info_handler(node n, GML_pair *read) {
  VERBOSESHOWER(8) << "Load node ... " << n;
  struct GML_pair *tmp = read;

  double x = 0;
  double y = 0;
  double z = 0;
  double t = 0;

  while (tmp) {
    //printf ("*KEY* : %s %f \n", tmp->key, tmp->value.floating);
    if (((string)(tmp->key)).find("x") < 1)
      x = tmp->value.floating;
    if (((string)(tmp->key)).find("y") < 1)
      y = tmp->value.floating;
    if (((string)(tmp->key)).find("z") < 1)
      z = tmp->value.floating;
    if (((string)(tmp->key)).find("t") < 1)
      t = tmp->value.floating;

    tmp = tmp->next;
  }

  vMap[n] = make_shared<Vertex>(x, y, z, t);
}

// use with graphviz (on Mac: brew install graphviz --with-app)
// dot GVfile.gv -Tpdf -o outputPDF.pdf

void PartonShower::SaveAsGV(string fName) {
  ofstream gv;
  gv.open(fName.c_str());

  // Simple directed graph left->right in dot/gv format for usage with graphviz ...
  // nodes show (time) and arrows (pT)
  gv << "digraph \"graph\" {" << endl;
  gv << endl;
  gv << "rankdir=\"LR\";" << endl;
  gv << "node [shape=plaintext, fontsize=11];"
     << endl; //, shape=circle]; //plaintext];
  gv << "edge [fontsize=10];" << endl;
  gv << endl;
  //gv<<"0 -> 1"<<endl;
  node_iterator nIt, nEnd;

  int n = 0;
  string label;
  string label2;

  for (nIt = nodes_begin(), nEnd = nodes_end(); nIt != nEnd; ++nIt) {
    label = ("[label=\"" + to_string(n) + "(");

    label2 = ")\"];";
    stringstream stream;

    stream << fixed << setprecision(2) << (vMap[*nIt]->x_in().t());
    gv << n << " " << label << stream.str() << label2 << endl;

    n++;
  }

  gv << endl;

  edge_iterator eIt, eEnd;

  for (eIt = edges_begin(), eEnd = edges_end(); eIt != eEnd; ++eIt) {
    //label = ("[label=\"(");
    if ((pMap[*eIt]->pstat()) == -13) {
      // missing energy-momentum
      label = ("[style=\"dotted\" color=\"red\"label=\"(");
    } else if ((pMap[*eIt]->pstat()) == -11) {
      // liquefied partons
      label = ("[color=\"red\"label=\"(");
    } else if ((pMap[*eIt]->pstat()) == -17) {
      // thermal partons draw from the medium (negative partons)
      label = ("[color=\"darkorange\"label=\"(");
    } else if ((pMap[*eIt]->pstat()) == -1) {
      // thermal partons draw from the medium (negative partons)
      label = ("[color=\"green\"label=\"(");
    } else if ((pMap[*eIt]->t()) < 4.0) {
      // small virtuality parton
      label = ("[color=\"blue\"label=\"(");
    } else {
      label = ("[label=\"(");
    }
    label2 = ")\"];";
    stringstream stream;
    if ((pMap[*eIt]->pstat()) == -13) {
      stream << std::scientific << setprecision(1) << (pMap[*eIt]->e()) << ","
             << (pMap[*eIt]->pt()) << "," << (pMap[*eIt]->t()) << ","
             << (pMap[*eIt]->pid());
    } else {
      stream << fixed << setprecision(2) << (pMap[*eIt]->e()) << ","
             << (pMap[*eIt]->pt()) << "," << (pMap[*eIt]->t()) << ","
             << (pMap[*eIt]->pid());
    }

    gv << (to_string(eIt->source().id()) + "->" + to_string(eIt->target().id()))
       << " " << label << stream.str() << label2 << endl;
  }
  gv << endl;
  gv << "}" << endl;
  gv.close();
}

void PartonShower::SaveAsGraphML(string fName) {
  // Think about using tinyxml2 in future (if needed ...)

  ofstream g;
  g.open(fName.c_str());

  g << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
  g << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"" << endl;
  g << "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"" << endl;
  g << "xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns" << endl;
  g << " http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">" << endl;

  g << "<key id=\"nlabel\" for=\"node\" attr.name=\"label\" "
       "attr.type=\"string\"/>"
    << endl;
  g << "<key id=\"nx\" for=\"node\" attr.name=\"x\" attr.type=\"double\"/>"
    << endl;
  g << "<key id=\"ny\" for=\"node\" attr.name=\"y\" attr.type=\"double\"/>"
    << endl;
  g << "<key id=\"nz\" for=\"node\" attr.name=\"z\" attr.type=\"double\"/>"
    << endl;
  g << "<key id=\"nt\" for=\"node\" attr.name=\"t\" attr.type=\"double\"/>"
    << endl;

  g << "<key id=\"elabel\" for=\"edge\" attr.name=\"label\" "
       "attr.type=\"string\"/>"
    << endl;
  g << "<key id=\"epl\" for=\"edge\" attr.name=\"plabel\" attr.type=\"int\"/>"
    << endl;
  g << "<key id=\"epid\" for=\"edge\" attr.name=\"pid\" attr.type=\"int\"/>"
    << endl;
  g << "<key id=\"estat\" for=\"edge\" attr.name=\"pstat\" attr.type=\"int\"/>"
    << endl;

  g << "<key id=\"ept\" for=\"edge\" attr.name=\"pT\" attr.type=\"double\"/>"
    << endl;
  g << "<key id=\"eeta\" for=\"edge\" attr.name=\"eta\" attr.type=\"double\"/>"
    << endl;
  g << "<key id=\"ephi\" for=\"edge\" attr.name=\"phi\" attr.type=\"double\"/>"
    << endl;
  g << "<key id=\"ee\" for=\"edge\" attr.name=\"E\" attr.type=\"double\"/>"
    << endl;

  g << "<graph id=\"G\" edgedefault=\"directed\">" << endl;

  node_iterator nIt, nEnd;
  int n = 0;

  for (nIt = nodes_begin(), nEnd = nodes_end(); nIt != nEnd; ++nIt) {
    stringstream stream;

    stream << fixed << setprecision(2) << (vMap[*nIt]->x_in().t());

    g << "<node id=\"" << n << "\">" << endl;
    //g<<"<data key=\"nlabel\">"<<to_string(n)+"("+to_string(vMap[*nIt]->x_in().t())+")"<<"</data>"<<endl;
    g << "<data key=\"nlabel\">" << to_string(n) + "(" + stream.str() + ")"
      << "</data>" << endl;
    g << "<data key=\"nx\">" << vMap[*nIt]->x_in().x() << "</data>" << endl;
    g << "<data key=\"ny\">" << vMap[*nIt]->x_in().y() << "</data>" << endl;
    g << "<data key=\"nz\">" << vMap[*nIt]->x_in().z() << "</data>" << endl;
    g << "<data key=\"nt\">" << vMap[*nIt]->x_in().t() << "</data>" << endl;
    g << "</node>" << endl;
    n++;
  }

  edge_iterator eIt, eEnd;
  n = 0;
  for (eIt = edges_begin(), eEnd = edges_end(); eIt != eEnd; ++eIt) {
    g << "<edge id=\"" << n << "\" source=\"" << to_string(eIt->source().id())
      << "\" target=\"" << to_string(eIt->target().id()) << "\">" << endl;
    g << "<data key=\"elabel\">" << to_string((pMap[*eIt]->pt())) << "</data>"
      << endl;
    g << "<data key=\"epl\">" << pMap[*eIt]->plabel() << "</data>" << endl;
    g << "<data key=\"epid\">" << pMap[*eIt]->pid() << "</data>" << endl;
    g << "<data key=\"estat\">" << pMap[*eIt]->pstat() << "</data>" << endl;
    g << "<data key=\"ept\">" << pMap[*eIt]->pt() << "</data>" << endl;
    g << "<data key=\"eeta\">" << pMap[*eIt]->eta() << "</data>" << endl;
    g << "<data key=\"ephi\">" << pMap[*eIt]->phi() << "</data>" << endl;
    g << "<data key=\"ee\">" << pMap[*eIt]->e() << "</data>" << endl;
    g << "</edge>" << endl;
    n++;
  }

  g << "</graph>" << endl;
  g << "</graphml>" << endl;

  g.close();

} // end namespace Jetscape
} // namespace Jetscape
