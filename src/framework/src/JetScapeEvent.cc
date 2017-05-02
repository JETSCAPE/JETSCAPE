// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#include "JetScapeEvent.h"
#include<iostream>

using namespace std;

JetScapeEvent::JetScapeEvent()
{
}

JetScapeEvent::JetScapeEvent(const JetScapeEvent &c){
    partonCollection.clear();
    const vector<Parton> tmp = c.getPartonCollection();
    for(unsigned int ipart=0; ipart<tmp.size(); ipart++){
        partonCollection.push_back(c.getParton(ipart));
    }
}

JetScapeEvent::~JetScapeEvent()
{
    partonCollection.clear();
}

const vector<Parton>& JetScapeEvent::getPartonCollection() const {
    return partonCollection;
}

const Parton& JetScapeEvent::getParton(int idx) const {
    return partonCollection.at(idx);
}

void JetScapeEvent::addParton(Parton &p){
    partonCollection.push_back(p);
}

void JetScapeEvent::addPartonShower(shared_ptr<PartonShower> ps){
    for(unsigned int ipart=0; ipart<ps->GetNumberOfPartons(); ipart++){
        partonCollection.push_back(*(ps->GetPartonAt(ipart)));
    } 
}

void JetScapeEvent::deleteParton(int idx){
    partonCollection.erase(partonCollection.begin()+idx); //inefficient delete!!
}
