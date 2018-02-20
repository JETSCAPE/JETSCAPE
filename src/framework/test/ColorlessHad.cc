#include "ColorlessHad.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>

using namespace Jetscape;
using namespace Pythia8;

Pythia pythia;
Event& event      = pythia.event;
ParticleData& pdt = pythia.particleData;

//Hadrons output file
ofstream hadfile;

ColorlessHad::ColorlessHad()
{
  SetId("ColorlessHad");
  VERBOSE(8);
}

ColorlessHad::~ColorlessHad()
{
  VERBOSE(8);
}

void ColorlessHad::Init()
{
  //Open output file
  hadfile.open("CH_myhad.dat");

  JSDEBUG<<"Initialize ColorlessHad";
  VERBOSE(8);

  // No event record printout.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Standard settings
  pythia.readString("ProcessLevel:all = off");
  
  // Don't let pi0 decay
  pythia.readString("111:mayDecay = off");

  // XML settings to be incorporated

  // And initialize
  pythia.init();

}

void ColorlessHad::WriteTask(weak_ptr<JetScapeWriter> w)
{
   VERBOSE(8);
   w.lock()->WriteComment("Hadronization Module : "+GetId());
}

void ColorlessHad::DoHadronization(vector<vector<shared_ptr<Parton>>>& shower, vector<shared_ptr<Hadron>>& hOut, vector<shared_ptr<Parton>>& pOut)
{
  INFO<<"Start Hadronizing using PYTHIA Lund string model (does NOT use color flow, needs to be tested)...";
  event.reset();

  //Hadronize each shower separately
  for(unsigned int ishower=0; ishower <  shower.size(); ++ishower)
  {
    vector<shared_ptr<Parton>>& pIn = shower.at(ishower);
    JSDEBUG<<"Shower#"<<ishower+1 << ". Number of partons to hadronize: " << pIn.size();

    int col[pIn.size()+2], acol[pIn.size()+2], isdone[pIn.size()+2];
    memset( col, 0, (pIn.size()+2)*sizeof(int) ), memset( acol, 0, (pIn.size()+2)*sizeof(int) ), memset( isdone, 0, (pIn.size()+2)*sizeof(int) );
  
    //Find number of quarks
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
  
    //Find number of strings
    int nstrings=max(int(double(nquarks)/2.+0.6),1);
    JSDEBUG << "#Strings = " << nstrings;
    //If there are no quarks, need to attach two of them
    int istring=0;
    int one_end[nstrings], two_end[nstrings];
    if (nquarks==0) {
      //First quark
      FourVector p1(0.,0.,1000.,1001.);
      FourVector x1;
      pIn.push_back(std::make_shared<Parton> (0,1,0,p1,x1));
      isquark[nquarks]=pIn.size()-1;
      nquarks+=1;
      isdone[pIn.size()-1]=1;
      one_end[0]=pIn.size()-1;
      INFO << "Attached quark remnant flying down +Pz beam";
      //Second quark
      FourVector p2(0.,0.,-1000.,1001.);
      FourVector x2;
      pIn.push_back(std::make_shared<Parton> (0,1,0,p2,x2));
      isquark[nquarks]=pIn.size()-1;
      nquarks+=1;
      isdone[pIn.size()-1]=1;
      two_end[istring]=pIn.size()-1;
      INFO << "Attached quark remnant flying down -Pz beam";
    }

    //Assign ends of strings (order matters in this algo)
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
          FourVector p(0.,0.,1000.,1001.);
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

    //Assign gluons to a certain string
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

    //Build up chain using gluons assigned to each string, in a closest pair order
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
      //Attach second end
      if (col[link]==lab_col-1) col[two_end[ns]]=0, acol[two_end[ns]]=lab_col-1;
      else col[two_end[ns]]=lab_col-1, acol[two_end[ns]]=0;
    }
    //Changing identity of quarks to be consistent with color charge
    for( int iq=0; iq <  nquarks; ++iq)
    {
      if (col[isquark[iq]]!=0) { 
        if (pIn[isquark[iq]]->pid()<0) pIn[isquark[iq]]->set_id(-pIn[isquark[iq]]->pid());
      }
      else {
        if (pIn[isquark[iq]]->pid()>0) pIn[isquark[iq]]->set_id(-pIn[isquark[iq]]->pid());
      }
    }

    //Introduce partons into PYTHIA
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
        //hOut.push_back(std::make_shared<Hadron> (Hadron (0,ide,0,p,x)));
        //INFO << "Produced Hadron has id = " << pythia.event[ipart].id();
        //Print on output file
        hadfile << pythia.event[ipart].px() << " " << pythia.event[ipart].py() << " " << pythia.event[ipart].pz() << " " << pythia.event[ipart].e() << " " << pythia.event[ipart].id() << " " << pythia.event[ipart].charge() << endl;
      }
    } 
    INFO<<"#Showers done: " << ishower+1 << ". There are " << hOut.size() << " hadrons and " << pOut.size() << " partons after PYTHIA Hadronization";
    hadfile << "NEXT" << endl;
  } // End shower loop
}
