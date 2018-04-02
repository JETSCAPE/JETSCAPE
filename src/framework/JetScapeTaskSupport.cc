/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * For the full list of contributors see AUTHORS.
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#include "JetScapeTaskSupport.h"
#include "JetScapeLogger.h"
#include <stdlib.h>

using namespace std;

namespace Jetscape {

  // static member initialization
  JetScapeTaskSupport* JetScapeTaskSupport::m_pInstance = nullptr;
  shared_ptr<std::mt19937> JetScapeTaskSupport::one_for_all_ = nullptr;
  unsigned int JetScapeTaskSupport::random_seed_=0;
  bool JetScapeTaskSupport::initialized_=false;
  bool JetScapeTaskSupport::one_generator_per_task_=false;
  
  // ---------------------------------------------------------------------------
  JetScapeTaskSupport* JetScapeTaskSupport::Instance()  {
    if (!m_pInstance) {
      m_pInstance = new JetScapeTaskSupport;
      VERBOSE(1)<<"Created JetScapeTaskSupport Instance";
    }    
    return m_pInstance;
  }
  
  // ---------------------------------------------------------------------------
  int JetScapeTaskSupport::RegisterTask(){
    VERBOSE(1) << "JetScapeTaskSupport::RegisterTask called, answering " << current_task_number_;
    current_task_number_++;
    return current_task_number_-1;
  }
  
  // ---------------------------------------------------------------------------
  void JetScapeTaskSupport::ReadSeedFromXML(){
    VERBOSE(1) << "JetScapeTaskSupport::ReadSeedFromXML called. ";

    // xml limits us to unsigned int :-/ -- but so does 32 bits Mersenne Twist
    tinyxml2::XMLElement *RandomXmlDescription=JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Random" );
    tinyxml2::XMLElement *xmle=0;
    if ( RandomXmlDescription ){
      xmle = RandomXmlDescription->FirstChildElement( "seed" ); if ( !xmle ) throw std::runtime_error("Cannot parse xml");
      xmle->QueryUnsignedText(&random_seed_);
    } else {
      WARN << "No <Random> element found in xml, seeding to 0";
    }
    
    VERBOSE(7) <<"Seeding JetScapeTaskSupport to "<< random_seed_ ;
    if ( random_seed_==0 ){
      // random seed
      random_seed_ = std::chrono::system_clock::now().time_since_epoch().count();
      VERBOSE(7) <<"JetScapeTaskSupport found seed 0, using one engine for all and reseeding to "<< random_seed_ ;
      one_generator_per_task_=false;
    } else {
      VERBOSE(7) <<"JetScapeTaskSupport found seed " << random_seed_ << ", using individual engines with seeds created from "<< random_seed_ ;
      one_generator_per_task_=true;
    }

    one_for_all_ = make_shared<std::mt19937>(random_seed_);
	  
    // VERBOSE(7) << "Setting random seed for mt19937 to " << seed;
    // generator.seed(seed);
    // ZeroOneDistribution = uniform_real_distribution<double> { 0.0, 1.0 };

    initialized_=true;      
  }

  // ---------------------------------------------------------------------------
  shared_ptr<std::mt19937> JetScapeTaskSupport::get_mt19937_generator( int TaskId ) {
    if (!initialized_){
      WARN << "Trying to use JetScapeTaskSupport::get_mt19937_generator before initialization";
      throw std::runtime_error("Trying to use JetScapeTaskSupport::get_mt19937_generator before initialization");
    }

    if ( !one_for_all_ ){
      throw std::runtime_error("generator not initialized?");
    }

    // In this case, we use our own engine to create a repeatable
    // sequence of seeds and hand over a properly seeded engine
    if ( one_generator_per_task_ ){
      if ( random_seed_==0 ) throw std::runtime_error("This should never happen");
      // reseed to be on the safe side
      one_for_all_->seed( random_seed_ );
      // Advance according to TaskId
      // Note that this method can be lied to.
      // Could design a safer interface but for now, trust the user
      one_for_all_->discard(TaskId);
      // And get the unique seed for this task
      unsigned int localseed = (*one_for_all_)();
      JSDEBUG << "Asked by " << TaskId << " for an individual generator, returning one seeded with " << localseed;
      return make_shared<std::mt19937>( localseed );
    }

    // this singleton owns the generator(s) and keeps them until deletion
    JSDEBUG << "Asked by " << TaskId << " for the static generator, returning one originally seeded with " << random_seed_;
    return one_for_all_;
  }
  
  // ---------------------------------------------------------------------------
  
} // end namespace Jetscape
