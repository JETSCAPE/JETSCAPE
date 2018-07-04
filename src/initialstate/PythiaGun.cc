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

// Create a pythia collision at a specified point and return the two inital hard partons

#include "PythiaGun.h"
#include <sstream>

using namespace std;


PythiaGun::~PythiaGun()
{
  VERBOSE(8);
}

void PythiaGun::InitTask()
{

  JSDEBUG<<"Initialize PythiaGun"; 
  VERBOSE(8);

  // Show initialization at INFO level
  readString("Init:showProcesses = off");
  readString("Init:showChangedSettings = off");
  readString("Init:showMultipartonInteractions = off");
  readString("Init:showChangedParticleData = off");
  if ( JetScapeLogger::Instance()->GetInfo() ) {
    readString("Init:showProcesses = on");
    readString("Init:showChangedSettings = on");
    readString("Init:showMultipartonInteractions = on");
    readString("Init:showChangedParticleData = on");
  }

  // No event record printout.
  readString("Next:numberShowInfo = 0"); 
  readString("Next:numberShowProcess = 0"); 
  readString("Next:numberShowEvent = 0"); 

  // Standard settings
<<<<<<< HEAD
  readString("HardQCD:all = off"); // will repeat this line in the xml for demonstration  
=======
  readString("HardQCD:all = on"); // will repeat this line in the xml for demonstration  
>>>>>>> a8fdc27b03dd460fc82996bb1aa469ebf9cbe306
  readString("HadronLevel:Decay = off");
  readString("HadronLevel:all = off");
  readString("PartonLevel:ISR = on");
  readString("PartonLevel:MPI = on");
<<<<<<< HEAD
  readString("PartonLevel:FSR = on");
=======
  readString("PartonLevel:FSR = off");
>>>>>>> a8fdc27b03dd460fc82996bb1aa469ebf9cbe306
  readString("PromptPhoton:all=on");
  readString("WeakSingleBoson:all=off");
  readString("WeakDoubleBoson:all=off");

  tinyxml2::XMLElement *PythiaXmlDescription=GetHardXML()->FirstChildElement("PythiaGun");
  tinyxml2::XMLElement *xmle; 

  string s;
  // For parsing text
  stringstream numbf(stringstream::app|stringstream::in|stringstream::out);
  numbf.setf(ios::fixed, ios::floatfield);       numbf.setf(ios::showpoint);       numbf.precision(1);
  stringstream numbi(stringstream::app|stringstream::in|stringstream::out);
    
  if ( !PythiaXmlDescription ) {
    WARN << "Cannot initialize Pythia Gun";
    throw std::runtime_error("Cannot initialize Pythia Gun");
  }

  xmle = PythiaXmlDescription->FirstChildElement( "name" ); if ( !xmle ) throw std::runtime_error("Cannot parse xml");
  s = xmle->GetText();
  SetId(s);
  // cout << s << endl;
  
  xmle = PythiaXmlDescription->FirstChildElement( "pTHatMin" ); if ( !xmle ) throw std::runtime_error("Cannot parse xml");
  xmle->QueryDoubleText(&pTHatMin);
  xmle = PythiaXmlDescription->FirstChildElement( "pTHatMax" ); if ( !xmle ) throw std::runtime_error("Cannot parse xml");
  xmle->QueryDoubleText(&pTHatMax);
  
  VERBOSE(7) <<"Pythia Gun with "<< pTHatMin << " < pTHat < " << pTHatMax ;
  
  numbf.str("PhaseSpace:pTHatMin = "); numbf << pTHatMin;
  readString ( numbf.str() );
  numbf.str("PhaseSpace:pTHatMax = "); numbf << pTHatMax;
  readString ( numbf.str() );
  
  // random seed
  // xml limits us to unsigned int :-/ -- but so does 32 bits Mersenne Twist
  tinyxml2::XMLElement *RandomXmlDescription=JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Random" );
  readString("Random:setSeed = on");
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
  readString( numbi.str() );
  
  // Species
  readString("Beams:idA = 2212");
  readString("Beams:idB = 2212");
  
  // Energy
  xmle = PythiaXmlDescription->FirstChildElement( "eCM" ); if ( !xmle ) throw std::runtime_error("Cannot parse xml");
  xmle->QueryDoubleText(&eCM);   
  numbf.str("Beams:eCM = "); numbf << eCM;
  readString ( numbf.str() );
  
  xmle = PythiaXmlDescription->FirstChildElement( "LinesToRead" );
  if ( xmle ){
    std::stringstream lines;
    lines << xmle->GetText();
    int i=0;
    while(std::getline(lines,s,'\n')){
      if( s.find_first_not_of (" \t\v\f\r") == s.npos ) continue; // skip empty lines
      VERBOSE(7) <<  "Also reading in: " << s;
      readString (s);
    }
  }
  
  // And initialize
  init(); // Pythia>8.1
  
}

void PythiaGun::Exec()
{
  INFO<<"Run Hard Process : "<<GetId()<< " ...";
  VERBOSE(8)<<"Current Event #"<<GetCurrentEvent();
  //Reading vir_factor from xml for MATTER
<<<<<<< HEAD
  tinyxml2::XMLElement *eloss= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" );
  
=======
   tinyxml2::XMLElement *eloss= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" );
>>>>>>> a8fdc27b03dd460fc82996bb1aa469ebf9cbe306
  if ( !eloss ) {
    WARN << "Couldn't find tag Eloss";
    throw std::runtime_error ("Couldn't find tag Eloss");    
  }
  tinyxml2::XMLElement *matter=eloss->FirstChildElement("Matter");
  if ( !matter ) {
    WARN << "Couldn't find tag Eloss -> Matter";
    throw std::runtime_error ("Couldn't find tag Eloss -> Matter");
  }
  double vir_factor;
  matter->FirstChildElement("vir_factor")->QueryDoubleText(&vir_factor);

  bool flag62=false;
  vector<Pythia8::Particle> p62;
<<<<<<< HEAD
  vector<int> hq_channels;
  vector<int> hq_mothers;
=======
  
>>>>>>> a8fdc27b03dd460fc82996bb1aa469ebf9cbe306

  // sort by pt
  struct greater_than_pt {
    inline bool operator() (const Pythia8::Particle& p1, const Pythia8::Particle& p2) {
      return ( p1.pT() > p2.pT());
    }
  };

  do{
    next();
    p62.clear();
<<<<<<< HEAD
    hq_channels.clear();
    hq_mothers.clear();
=======
>>>>>>> a8fdc27b03dd460fc82996bb1aa469ebf9cbe306
    
    // pTarr[0]=0.0; pTarr[1]=0.0;
    // pindexarr[0]=0; pindexarr[1]=0;

    for(int parid=0; parid<event.size(); parid++){
      if ( parid<3 )continue;      // 0, 1, 2: total event and beams      
      Pythia8::Particle& particle = event[parid];
<<<<<<< HEAD
 
      //INFO<<"id: "<<particle.id()<<" , status: "<<particle.status();

      // only accept particles after MPI
      if ( particle.status()!=62 ) continue;
      
      // only accept gluons and quarks
      // Also accept Gammas to put into the hadron's list
      // put heavy quarks in!!!

      if ( fabs( particle.id() ) > 5 && (particle.id() !=21 && particle.id() !=22) ) continue;

=======

      // only accept particles after MPI
      if ( particle.status()!=62 ) continue;
      // only accept gluons and quarks
      // Also accept Gammas to put into the hadron's list
      if ( fabs( particle.id() ) > 3 && (particle.id() !=21 && particle.id() !=22) ) continue;
      
>>>>>>> a8fdc27b03dd460fc82996bb1aa469ebf9cbe306
      // reject rare cases of very soft particles that don't have enough e to get
      // reasonable virtuality
      if ( particle.pT() < 1.0/sqrt(vir_factor) ) continue;
 	
	//if(particle.id()==22) cout<<"########this is a photon!######" <<endl;
      // accept
      p62.push_back( particle );
<<<<<<< HEAD
      
      //specify the production type of heavy quark
      if (abs(particle.id()) == 4 || abs(particle.id()) == 5)
	    {
                //INFO << "particle: " << particle.id() << " event_id: " << parid << " status: " <<particle.status();
		            int copy=particle.iTopCopyId();
                int mother1=event[copy].mother1();
                int mother2=event[copy].mother2();
                if(mother1>0 && mother2>0)
		  {
 			INFO << "generated by flavor creation";
                        hq_channels.push_back(1);
		  }
                else
                {
		        int daughter1=event[mother1].daughter1();
                        int daughter2=event[mother1].daughter2();
                        if(abs(event[daughter1].id())!=particle.id()||abs(event[daughter2].id())!=particle.id())
						          	{
								          hq_channels.push_back(-1);
								          continue;
							          }
                        bool has_scat=false;
                                                        while(event[daughter1].status()<0)
							{
								int tempdaughter1 = event[daughter1].daughter1();
								int tempdaughter2 = event[daughter1].daughter2();
								if (event[daughter1].status() == -31 || event[daughter1].status() == -21)
								{
									has_scat=true;
									break;
								}
								
								if(event[tempdaughter1].id()==event[daughter1].id())
								{
									daughter1=tempdaughter1;
								}
								else
								{
									daughter1=tempdaughter2;
								}
							}
                                                        while(event[daughter2].status()<0)
							{
								int tempdaughter1 = event[daughter2].daughter1();
								int tempdaughter2 = event[daughter2].daughter2();
								if (event[daughter2].status() == -31 || event[daughter2].status() == -21)
								{
									has_scat=true;
									break;

								}
								if(event[tempdaughter1].id()==event[daughter2].id())
								{
									daughter2=tempdaughter1;
								}
								else
								{
									daughter2=tempdaughter2;
								}
							}
      if(has_scat==true)
			{
				INFO << "generated by flavor excitation";
   				hq_channels.push_back(2);
			}
      else
			{
				INFO << "generated by gluon splitting";
				hq_channels.push_back(3);
			}

	  }
		hq_mothers.push_back(mother1);
	}
      else
        {
    hq_channels.push_back(-1);
		hq_mothers.push_back(particle.mother1());
        }
=======
>>>>>>> a8fdc27b03dd460fc82996bb1aa469ebf9cbe306

    }

    // if you want at least 2
    if ( p62.size() < 2 ) continue;
    //if ( p62.size() < 1 ) continue;
    
    // Now have all candidates, sort them
    // sort by pt
<<<<<<< HEAD
    //std::sort( p62.begin(), p62.end(), greater_than_pt() );
=======
    std::sort( p62.begin(), p62.end(), greater_than_pt() );
>>>>>>> a8fdc27b03dd460fc82996bb1aa469ebf9cbe306
    // // check...
    // for (auto& p : p62 ) cout << p.pT() << endl;
    
    flag62=true;

  }while(!flag62);


  double p[4], xLoc[4];

  // This location should come from an initial state
  for (int i=0;i<=3; i++) {
    xLoc[i] = 0.0;
  };

  // // Roll for a starting point
  // // See: https://stackoverflow.com/questions/15039688/random-generator-from-vector-with-probability-distribution-in-c
  // std::random_device device;
  // std::mt19937 engine(device()); // Seed the random number engine

  
  if (!ini) {
      WARN << "No initial state module, setting the starting location to 0. Make sure to add e.g. trento before PythiaGun.";
  } else {
    auto num_bin_coll = ini->GetNumOfBinaryCollisions();
    if ( num_bin_coll.size()==0 ){
      WARN << "num_of_binary_collisions is empty, setting the starting location to 0. Make sure to add e.g. trento before PythiaGun.";
    } else {	 
      std::discrete_distribution<> dist( begin(num_bin_coll),end(num_bin_coll) ); // Create the distribution
    
      // Now generate values
      auto idx = dist( *GetMt19937Generator() );
      auto coord = ini->CoordFromIdx( idx );
      xLoc[1] = get<0>( coord );
      xLoc[2] = get<1>( coord );
    }
  }
    

  // Loop through particles

  // Only top two
  //for(int np = 0; np<2; ++np){

  // Accept them all

  int hCounter = 0 ;
  for(int np = 0; np<p62.size(); ++np){
    Pythia8::Particle& particle = p62.at( np );

<<<<<<< HEAD
    //VERBOSE(7)
    VERBOSE(7) <<"Adding particle with pid = " << particle.id()
=======
    VERBOSE(7)<<"Adding particle with pid = " << particle.id()
>>>>>>> a8fdc27b03dd460fc82996bb1aa469ebf9cbe306
	      <<" at x=" << xLoc[1]
	      <<", y=" << xLoc[2]
	      <<", z=" << xLoc[3];
    
    VERBOSE(7) <<"Adding particle with pid = " << particle.id()
	       << ", pT = " << particle.pT()
	       << ", y = " << particle.y()
	       << ", phi = " << particle.phi()
	       << ", e = " << particle.e();

    VERBOSE(7) <<" at x=" << xLoc[1]
	       <<", y=" << xLoc[2]
	       <<", z=" << xLoc[3];
    if(particle.id() !=22)
    {
<<<<<<< HEAD
        AddParton(make_shared<Parton>(0, particle.id(),0,particle.pT(),particle.y(),particle.phi(),particle.e(),xLoc));
        //special treatment for heavy quarks
        if (abs(particle.id()) == 4 || abs(particle.id()) == 5) 
	{
                shared_ptr<Parton> hq=GetPartonAt(GetNHardPartons()-1);
       		      hq->set_hq_channel(hq_channels[np]);
                hq->set_hq_mother_id(hq_mothers[np]);
	}
=======
        AddParton(make_shared<Parton>(0, particle.id(),0,particle.pT(),particle.y(),particle.phi(),particle.e(),xLoc) );
>>>>>>> a8fdc27b03dd460fc82996bb1aa469ebf9cbe306
    }
    else
    {
	      AddHadron(make_shared<Hadron>(hCounter,particle.id(),particle.status(),particle.pT(),particle.eta(),particle.phi(),particle.e(),xLoc));
	      hCounter++;
    }
  }
  

 
  VERBOSE(8)<<GetNHardPartons();
}

