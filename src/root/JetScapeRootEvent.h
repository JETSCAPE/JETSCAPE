// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#ifndef ROOT_JetScapeRootEvent
#define ROOT_JetScapeRootEvent

#include "JetScapeRootPartonShower.h"
#include "TObject.h"

using namespace std;
using namespace Jetscape;

class JetScapeRootEvent : public TObject
  {
    
  public:
    
    JetScapeRootEvent();    
    ~JetScapeRootEvent();
    
    void SetEventNumber(int m_eventNumber) {eventNumber=m_eventNumber;}

    void AddPartonShower(unique_ptr<JetScapeRootPartonShower> mPS) {psVec.push_back(move(mPS));}

    void AddHydroHistory() {};
    void AddHadrons() {};
    // more can be added here ...

    int GetEventNumber() {return eventNumber;}
    int GetNumberOfPartonShowers() {return psVec.size();}
 
    unique_ptr<JetScapeRootPartonShower> GetPartonShowerAt(int i) {return move(psVec[i]);}       
    
  private:
    
    int eventNumber;

    vector<unique_ptr<JetScapeRootPartonShower>> psVec;
    
    ClassDef(JetScapeRootEvent,1)
      
  };

#endif
