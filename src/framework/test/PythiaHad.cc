#include "PythiaHad.h"
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

PythiaHad::PythiaHad()
{
  SetId("PythiaHad");
  VERBOSE(8);
}

PythiaHad::~PythiaHad()
{
  VERBOSE(8);
}

void PythiaHad::Init()
{
  DEBUG<<"Initialize PythiaHad";
  VERBOSE(8);

  // No event record printout.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Standard settings
  pythia.readString("ProcessLevel:all = off");
  

  // ************
  // CANNOT MAKE THE XML WORK, GIVES SEG FAULT WHEN EXECUTING...
  // ************
  /*

  // XML settings
  tinyxml2::XMLElement *PythiaXmlDescription=GetHadronXML()->FirstChildElement("PythiaHad");
  tinyxml2::XMLElement *xmle;  

  string s;
  // For parsing text
  stringstream numbf(stringstream::app|stringstream::in|stringstream::out);
  numbf.setf(ios::fixed, ios::floatfield);       numbf.setf(ios::showpoint);       numbf.precision(1);
  stringstream numbi(stringstream::app|stringstream::in|stringstream::out);

  if ( !PythiaXmlDescription ) {
    WARN << "Cannot initialize Pythia Had";
    throw std::runtime_error("Cannot initialize Pythia Had");
  }

  xmle = PythiaXmlDescription->FirstChildElement( "name" ); if ( !xmle ) throw std::runtime_error("Cannot parse xml");
  s = xmle->GetText();
  SetId(s);
  cout << s << endl;

  //To use in future when choosing lund string parameters
  //xmle = PythiaXmlDescription->FirstChildElement( "pTHatMin" ); if ( !xmle ) throw std::runtime_error("Cannot parse xml");
  //xmle->QueryDoubleText(&pTHatMin);
  //xmle = PythiaXmlDescription->FirstChildElement( "pTHatMax" ); if ( !xmle ) throw std::runtime_error("Cannot parse xml");
  //xmle->QueryDoubleText(&pTHatMax);
  
  //VERBOSE(7) <<"Pythia Gun with "<< pTHatMin << " < pTHat < " << pTHatMax ;
  
  //numbf.str("PhaseSpace:pTHatMin = "); numbf << pTHatMin;
  //readString ( numbf.str() );
  //numbf.str("PhaseSpace:pTHatMax = "); numbf << pTHatMax;
  //readString ( numbf.str() );

  // random seed
  // xml limits us to unsigned int :-/ -- but so does 32 bits Mersenne Twist
  tinyxml2::XMLElement *RandomXmlDescription=JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Random" );
  pythia.readString("Random:setSeed = on");
  numbi.str("Random:seed = ");
  unsigned int seed = 0;
  if ( RandomXmlDescription ){
    xmle = RandomXmlDescription->FirstChildElement( "seed" ); if ( !xmle ) throw std::runtime_error("Cannot parse xml");
    xmle->QueryUnsignedText(&seed);
  } else {
    WARN << "No <Random> element found in xml, seeding to 0";
  }
  VERBOSE(7) <<"Seeding pythia to "<< seed ;
  numbi << seed;
  pythia.readString( numbi.str() );

  xmle = PythiaXmlDescription->FirstChildElement( "LinesToRead" );
  if ( xmle ){
    std::stringstream lines; lines << xmle->GetText();
    int i=0;
    while(std::getline(lines,s,'\n')){
      if( s.find_first_not_of (" \t\v\f\r") == s.npos ) continue; // skip empty lines
      VERBOSE(7) <<  "Also reading in: " << s;
      pythia.readString (s);
    }
  }
  
  */

  // And initialize
  pythia.init();

}

void PythiaHad::WriteTask(weak_ptr<JetScapeWriter> w)
{
   VERBOSE(8);
   w.lock()->WriteComment("Hadronization Module : "+GetId());
   w.lock()->WriteComment("Hadronization to be implemented accordingly ...");
}

void PythiaHad::DoHadronization(vector<shared_ptr<Parton>>& pIn, vector<shared_ptr<Hadron>>& hOut, vector<shared_ptr<Parton>>& pOut)
{
  INFO<<"Start Hadronizing using Lund string model...";
  event.reset();
  INFO<<"Number of partons to hadronize: " << pIn.size();

  int col[pIn.size()+1], acol[pIn.size()+1], isdone[pIn.size()+1];
  memset( col, 0, (pIn.size()+1)*sizeof(int) ), memset( acol, 0, (pIn.size()+1)*sizeof(int) ), memset( isdone, 0, (pIn.size()+1)*sizeof(int) );
  
  //Find number of quarks
  int nquarks=0;
  int isquark[pIn.size()+1];
  memset( isquark, 0, (pIn.size()+1)*sizeof(int) );
  for(unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
  {
    if (abs(pIn[ipart]->pid())<=6) {
      isquark[nquarks]=ipart;
      nquarks+=1;  
    }
  }
  INFO << "#Quarks = " << nquarks;
  
  //Find number of strings
  int nstrings=int(double(nquarks)/2.+0.6);
  INFO << "#Strings = " << nstrings;

  //Find ends of strings (order matters in this algo)
  int istring=0;
  int one_end[nstrings], two_end[nstrings];
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
        INFO << " Attached remnant flying down +Pz beam";
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
        cout << " Linked parton= " << next_link << endl;
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

  /*
  //Find closest partners:
  for(unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
  {
    //If it is a gluon, need to find two closest partners, if it is a quark only the closest
    int thisid=0;
    if (abs(pIn[ipart]->pid())<=6) thisid=1;
    else if (pIn[ipart]->pid()==21) thisid=2;
    else continue;
    double min_delR[2]={100000.};
    for(unsigned int jpart=0; jpart <  pIn.size(); ++jpart)
    {
      if (ipart == jpart) continue;
      fjcore::PseudoJet pf(pIn[jpart]->px(),pIn[jpart]->py(),pIn[jpart]->pz(),pIn[jpart]->e());  
      double delR = pIn[ipart]->delta_R(pf);
      if (delR < min_delR[0]) min_delR[1]=min_delR[0], min_delR[0]=delR, partner[ipart][1]=partner[ipart][0], partner[ipart][0]=jpart;
      else if (delR < min_delR[1]) min_delR[1]=delR, partner[ipart][1]=jpart;
      cout << " jpart= " << jpart << endl;
    }  
  }
  */

  for (unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
  {
    cout << "Parton #" << ipart << " is a " << pIn[ipart]->pid() << "with energy = " << pIn[ipart]->e() << " with phi= " << pIn[ipart]->phi() << " and has col= " << col[ipart] << " and acol= " << acol[ipart] << endl;
  }
  for(unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
  {  
    int ide=pIn[ipart]->pid();
    double px=pIn[ipart]->px();
    double py=pIn[ipart]->py();
    double pz=pIn[ipart]->pz();
    double ee=pIn[ipart]->e();
    double mm=pdt.m0(int(ide));
    ee=std::sqrt(px*px+py*py+pz*pz+mm*mm);
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
      hOut.push_back(std::make_shared<Hadron> (Hadron (0,ide,0,p,x)));
      INFO << "Produced Hadron has id = " << pythia.event[ipart].id();
    }
  } 
  INFO<<"There are " << hOut.size() << " Hadrons and " << pOut.size() << " partons after Hadronization";

}
