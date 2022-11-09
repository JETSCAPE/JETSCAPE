// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#ifndef ROOT_JetScapeRootVertex
#define ROOT_JetScapeRootVertex

#include "JetClass.h"
#include "TObject.h"

using namespace std;
using namespace Jetscape;

class JetScapeRootVertex : public TObject
  {
    
  public:
    
    JetScapeRootVertex(); 
    JetScapeRootVertex(Vertex& srp);
    ~JetScapeRootVertex();
    
    void SetNodeId(int m_NodeId) {id=m_NodeId;}
    
    int GetNodeId() {return id;}
    
  private:
    
    int id;
    double x,y,z,t;
    
    ClassDef(JetScapeRootVertex,1)
      
  };

#endif
