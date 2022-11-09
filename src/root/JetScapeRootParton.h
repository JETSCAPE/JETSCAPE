// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#ifndef ROOT_JetScapeRootParton
#define ROOT_JetScapeRootParton

#include "JetScapeParticles.h"
#include "TObject.h"

using namespace std;
using namespace Jetscape;

class JetScapeRootParton : public TObject //, public Parton
  {
    
  public:
    
    JetScapeRootParton(); // {}; 
    JetScapeRootParton(Parton& srp); // : Parton(srp) {};
    ~JetScapeRootParton();
    
    void SetSourceNode(int m_inNode) {sourceNode=m_inNode;}
    void SetTargetNode(int m_outNode) {targetNode=m_outNode;}
    
    int GetSourceNode() {return sourceNode;}
    int GetTargetNode() {return targetNode;}
    
  private:
    
    int targetNode, sourceNode;
    
    int label, pid, stat;
    double pt, eta, phi, e;
    
    ClassDef(JetScapeRootParton,1)
      
  };

#endif
