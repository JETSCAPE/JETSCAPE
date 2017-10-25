#ifndef HADRONIZATION_H
#define HADRONIZATION_H

#include "JetScapeModuleBase.h"
#include "JetClass.hpp"
#include "JetScapeWriter.h"
#include "PartonShower.h"
//#include "helper.h"
#include <vector>
#include <random>
#include "tinyxml2.h"
#include "JetScapeXML.h"

namespace Jetscape {

class Hadronization : public JetScapeModuleBase, public std::enable_shared_from_this<Hadronization>
{
public:

  Hadronization();
  virtual ~Hadronization();
  virtual shared_ptr<Hadronization> Clone() const {return nullptr;}
  virtual void Init();
  virtual void Exec();
  virtual void DoHadronization(vector<shared_ptr<Parton>>& pIn, vector<shared_ptr<Hadron>>& hOut, vector<shared_ptr<Parton>>& pOut){};
  virtual void WriteTask(weak_ptr<JetScapeWriter> w); 
  virtual void Clear();

  tinyxml2::XMLElement* GetHadronXML() {return fd;}

  sigslot::signal3<vector<shared_ptr<Parton>>&, vector<shared_ptr<Hadron>>&, vector<shared_ptr<Parton>>&, multi_threaded_local> TransformPartons;

  vector<shared_ptr<Hadron>> GetHadrons(){return outHadrons;}
  vector<shared_ptr<Parton>> GetOutPartons(){return outPartons;}
  
  void AddInPartons(vector<shared_ptr<Parton>> ip) {inPartons = ip;}

  void SetTransformPartonsConnected(bool m_TransformPartonsConnected) {TransformPartonsConnected=m_TransformPartonsConnected;}
  const bool GetTransformPartonsConnected() {return TransformPartonsConnected;}

private:

  vector<shared_ptr<Parton>> inPartons;
  vector<shared_ptr<Hadron>> outHadrons;
  vector<shared_ptr<Parton>> outPartons;
  void DoHadronize();

  bool TransformPartonsConnected;

  tinyxml2::XMLElement *fd;

};


}

#endif
