/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * Modular, task-based framework
 * Intial Design: Joern Putschke, Kolja Kauder (Wayne State University)
 * For the full list of contributors see AUTHORS.

 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#include "ColorlessHadronization.h"
#include "JetScapeXML.h"
#include "JetScapeLogger.h"
#include "tinyxml2.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>

using namespace Jetscape;
using namespace Pythia8;

// Hadrons output file
//ofstream hadfile;

// Initialize static helper here
Pythia8::Pythia ColorlessHadronization::pythia ("IntentionallyEmpty",false);


ColorlessHadronization::ColorlessHadronization()
{
  SetId("ColorlessHadronization");
  VERBOSE(8);
}

ColorlessHadronization::~ColorlessHadronization()
{
  VERBOSE(8);
}

void ColorlessHadronization::Init()
{
  // Open output file
  //hadfile.open("CH_myhad.dat");

  tinyxml2::XMLElement *hadronization= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("JetHadronization" );

  if ( !hadronization ) {
    WARN << "Couldn't find tag Jet Hadronization";
    throw std::runtime_error ("Couldn't find tag Jet Hadronization");
  }
  if (hadronization) {
    string s = hadronization->FirstChildElement( "name" )->GetText();
    JSDEBUG << s << " to be initializied ...";

    // Read sqrts to know remnants energies
    double p_read_xml = 10000 ;
    hadronization->FirstChildElement("eCMforHadronization")->QueryDoubleText(&p_read_xml);
    p_fake = p_read_xml;

    JSDEBUG<<"Initialize ColorlessHadronization";
    VERBOSE(8);

    // No event record printout.
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");

    // Standard settings
    pythia.readString("ProcessLevel:all = off");
  
    // Don't let pi0 decay
    pythia.readString("111:mayDecay = off");

    // Don't let any hadron decay
    //pythia.readString("HadronLevel:Decay = off");

    // And initialize
    pythia.init();
  }
}

void ColorlessHadronization::WriteTask(weak_ptr<JetScapeWriter> w)
{
  VERBOSE(8);
  auto f = w.lock();
  if ( !f ) return; 
  f->WriteComment("Hadronization Module : "+GetId());
}

void ColorlessHadronization::DoHadronization(vector<vector<shared_ptr<Parton>>>& shower, vector<shared_ptr<Hadron>>& hOut, vector<shared_ptr<Parton>>& pOut)
{
  INFO<<"Start Hadronizing using PYTHIA Lund string model (does NOT use color flow, needs to be tested)...";
  Event& event      = pythia.event;
  ParticleData& pdt = pythia.particleData;

  // Hadronize positive (status = 0) and negative (status = -1) partons in a different space 
  for (int want_pos=1; want_pos>=0; --want_pos)
  {
    event.reset();

    // Set remnants momentum
    double rempx=0.2;
    double rempy=0.2;
    double rempz=p_fake;
    double reme=std::sqrt(std::pow(rempx,2.)+std::pow(rempy,2.)+std::pow(rempz,2.));

    // Hadronize all showers together
    vector<shared_ptr<Parton>> pIn;
    for(unsigned int ishower=0; ishower <  shower.size(); ++ishower)
    {
      for(unsigned int ipart=0; ipart <  shower.at(ishower).size(); ++ipart)
      {
        if (shower.at(ishower).at(ipart)->pstat()==0 && want_pos==1) pIn.push_back(shower.at(ishower).at(ipart));  // Positive
        if (shower.at(ishower).at(ipart)->pstat()==-1 && want_pos==0) pIn.push_back(shower.at(ishower).at(ipart));  // Negative
      }
      JSDEBUG<<"Shower#"<<ishower+1 << ". Number of partons to hadronize so far: " << pIn.size();
    }
    if (want_pos==1) INFO<<"# Positive Partons to hadronize: " << pIn.size();
    else INFO<<"# Negative Partons to hadronize: " << pIn.size(); 

    // Check whether event is empty (specially important for negative partons case)
    if (pIn.size()==0) continue;

    int col[pIn.size()+2], acol[pIn.size()+2], isdone[pIn.size()+2];
    memset( col, 0, (pIn.size()+2)*sizeof(int) ), memset( acol, 0, (pIn.size()+2)*sizeof(int) ), memset( isdone, 0, (pIn.size()+2)*sizeof(int) );
  
    // Find number of quarks
    int nquarks=0;
    int isquark[pIn.size()+2];
    memset( isquark, 0, (pIn.size()+2)*sizeof(int) );
    for(unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
    {
      if (abs(pIn[ipart]->pid())<=6) {
        isquark[nquarks]=ipart;
        nquarks+=1;  
      }
    }
    JSDEBUG << "#Quarks = " << nquarks;
  
    // Find number of strings
    int nstrings=max(int(double(nquarks)/2.+0.6),1);
    JSDEBUG << "#Strings = " << nstrings;
    
    // If there are no quarks, need to attach two of them
    int istring=0;
    int one_end[nstrings], two_end[nstrings];
    if (nquarks==0) { // Only attach remnants if event is not empty
      // First quark
      FourVector p1(rempx,rempy,rempz,reme);
      FourVector x1;
      pIn.push_back(std::make_shared<Parton> (0,1,0,p1,x1));
      isquark[nquarks]=pIn.size()-1;
      nquarks+=1;
      isdone[pIn.size()-1]=1;
      one_end[0]=pIn.size()-1;
      INFO << "Attached quark remnant flying down +Pz beam";
      // Second quark
      FourVector p2(rempx,rempy,-rempz,reme);
      FourVector x2;
      pIn.push_back(std::make_shared<Parton> (0,1,0,p2,x2));
      isquark[nquarks]=pIn.size()-1;
      nquarks+=1;
      isdone[pIn.size()-1]=1;
      two_end[istring]=pIn.size()-1;
      INFO << "Attached quark remnant flying down -Pz beam";
    }

    // Assign ends of strings (order matters in this algo)
    for(unsigned int iquark=0; iquark<nquarks; iquark++) {
      if (isdone[isquark[iquark]]==0) {
        isdone[isquark[iquark]]=1;
        one_end[istring]=isquark[iquark];
        double min_delR=10000.;
        int partner=-2;
        for(unsigned int jquark=0; jquark<nquarks; jquark++) {  
          if (iquark==jquark) continue;
          int d_jquark=isquark[jquark];
          if (isdone[d_jquark]==0) {
            fjcore::PseudoJet pf(pIn[d_jquark]->px(),pIn[d_jquark]->py(),pIn[d_jquark]->pz(),pIn[d_jquark]->e());
            double delR = pIn[isquark[iquark]]->delta_R(pf);    
            if (delR<min_delR) min_delR=delR, partner=jquark;
          }
        }
        if (partner!=-2) {
          isdone[isquark[partner]]=1;
          two_end[istring]=isquark[partner];
          istring+=1;
        }
        else {
          FourVector p(rempx,rempy,rempz,reme);
          FourVector x;
          pIn.push_back(std::make_shared<Parton> (0,1,0,p,x));
          isquark[nquarks]=pIn.size()-1;
          nquarks+=1;
          isdone[pIn.size()-1]=1;
          two_end[istring]=pIn.size()-1;
          INFO << "Attached quark remnant flying down +Pz beam";
        }
      }
    }

    // Assign gluons to a certain string
    int my_string[pIn.size()];
    memset( my_string, 0, pIn.size()*sizeof(int) );
    for(unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
    {
      if (pIn[ipart]->pid()==21) {
        double min_delR=100000.;
        for (int ns=0; ns<nstrings; ns++)
        {         
          int fq=one_end[ns];
          int sq=two_end[ns];
          fjcore::PseudoJet pfq(pIn[fq]->px(),pIn[fq]->py(),pIn[fq]->pz(),pIn[fq]->e());
          double f_delR = pIn[ipart]->delta_R(pfq);
          fjcore::PseudoJet psq(pIn[sq]->px(),pIn[sq]->py(),pIn[sq]->pz(),pIn[sq]->e());
          double s_delR = pIn[ipart]->delta_R(psq);
          double delR=(f_delR+s_delR)/2.;
          if (delR<min_delR) my_string[ipart]=ns, min_delR=delR;
        }
      }
    }

    // Build up chain using gluons assigned to each string, in a closest pair order
    int lab_col=102;
    for (int ns=0; ns<nstrings; ns++)
    {
      int tquark=one_end[ns];
      if (pIn[tquark]->pid()>0) col[tquark]=lab_col;
      else acol[tquark]=lab_col;
      lab_col+=1;
      int link=tquark;
      int changes=1;
      do {
        changes=0;
        double min_delR=100000.;
        int next_link=0;
        for(unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
        {
          if (pIn[ipart]->pid()==21 && isdone[ipart]==0 && my_string[ipart]==ns)
          {
            changes=1;
            fjcore::PseudoJet pf(pIn[ipart]->px(),pIn[ipart]->py(),pIn[ipart]->pz(),pIn[ipart]->e());
            double delR = pIn[link]->delta_R(pf);
            if (delR<min_delR) min_delR=delR, next_link=ipart;
          }
        }
        if (changes==1)
        {
          isdone[next_link]=1;
          if (col[link]==lab_col-1) col[next_link]=lab_col, acol[next_link]=lab_col-1;
          else col[next_link]=lab_col-1, acol[next_link]=lab_col;
          lab_col+=1;
          JSDEBUG << " Linked parton= " << next_link;
          link=next_link;
        }
      } while (changes==1);
      // Attach second end
      if (col[link]==lab_col-1) col[two_end[ns]]=0, acol[two_end[ns]]=lab_col-1;
      else col[two_end[ns]]=lab_col-1, acol[two_end[ns]]=0;
    }
    // Changing identity of quarks to be consistent with color charge
    for( int iq=0; iq <  nquarks; ++iq)
    {
      if (col[isquark[iq]]!=0) { 
        if (pIn[isquark[iq]]->pid()<0) pIn[isquark[iq]]->set_id(-pIn[isquark[iq]]->pid());
      }
      else {
        if (pIn[isquark[iq]]->pid()>0) pIn[isquark[iq]]->set_id(-pIn[isquark[iq]]->pid());
      }
    }

    // Introduce partons into PYTHIA
    /*
    for (unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
    {
      JSDEBUG << "Parton #" << ipart << " is a " << pIn[ipart]->pid() << "with energy = " << pIn[ipart]->e() << " with phi= " << pIn[ipart]->phi() << " and has col= " << col[ipart] << " and acol= " << acol[ipart];
    }
    */
    for(unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
    {  
      int ide=pIn[ipart]->pid();
      double px=pIn[ipart]->px();
      double py=pIn[ipart]->py();
      double pz=pIn[ipart]->pz();
      double ee=pIn[ipart]->e();
      double mm=pdt.m0(int(ide));
      ee=std::sqrt(px*px+py*py+pz*pz+mm*mm);
      if (col[ipart]==0 && acol[ipart]==0 && (ide==21 || abs(ide)<=6)) {
        INFO<<"Stopping because of colorless parton trying to be introduced in PYTHIA string";
        exit(0);
      }
      event.append(int(ide),23,col[ipart],acol[ipart],px,py,pz,ee,mm);
    }
    pythia.next();
    for (unsigned int ipart=0; ipart < event.size(); ++ipart)
    {
      if (event[ipart].isFinal())
      {
        int ide=pythia.event[ipart].id();
        FourVector p(pythia.event[ipart].px(),pythia.event[ipart].py(),pythia.event[ipart].pz(),pythia.event[ipart].e());
        FourVector x;
        if (want_pos==1) hOut.push_back(std::make_shared<Hadron> (Hadron (0,ide,0,p,x)));  // Positive
        else hOut.push_back(std::make_shared<Hadron> (Hadron (0,ide,-1,p,x)));  // Negative
        //INFO << "Produced Hadron has id = " << pythia.event[ipart].id();
        // Print on output file
        //hadfile << pythia.event[ipart].px() << " " << pythia.event[ipart].py() << " " << pythia.event[ipart].pz() << " " << pythia.event[ipart].e() << " " << pythia.event[ipart].id() << " " << pythia.event[ipart].charge() << endl;
      }
    } 
    INFO<<"#Showers hadronized together: " << shower.size() << ". There are " << hOut.size() << " hadrons and " << pOut.size() << " partons after PYTHIA Hadronization";
    //hadfile << "NEXT" << endl;
  
  } // End of positive or negative loop
}
