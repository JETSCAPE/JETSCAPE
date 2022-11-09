// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#ifndef ROOT_JetScapeRootPartonShower
#define ROOT_JetScapeRootPartonShower

#include "JetScapeRootParton.h"
#include "JetScapeRootVertex.h"

#include "TObject.h"

using namespace std;
using namespace Jetscape;

class JetScapeRootPartonShower : public TObject //, public Parton
  {
    
  public:
    
    JetScapeRootPartonShower(); // {};    
    ~JetScapeRootPartonShower();

    void AddVertex(unique_ptr<JetScapeRootVertex> mV) {vVec.push_back(move(mV));}
    void AddParton(unique_ptr<JetScapeRootParton> mP) {pVec.push_back(move(mP));}

    int GetNumberOfVertices() {return vVec.size();}
    int GetNumberOfPartons() {return pVec.size();}

    unique_ptr<JetScapeRootVertex> GetVertexAt(int i) {return move(vVec[i]);}
    unique_ptr<JetScapeRootParton> GetPartonAt(int i) {return move(pVec[i]);}
    
  private:

    vector<unique_ptr<JetScapeRootParton>> pVec;
    vector<unique_ptr<JetScapeRootVertex>> vVec;
    
    ClassDef(JetScapeRootPartonShower,1)
      
  };

#endif
