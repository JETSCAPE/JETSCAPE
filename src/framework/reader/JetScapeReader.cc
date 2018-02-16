// -----------------------------------------
// // JetScape (modular/task) based framework
// // Intial Design: Joern Putschke (2017)
// //                (Wayne State University)
// // -----------------------------------------
// // License and Doxygen-like Documentation to be added ...
//

#include "JetScapeReader.h"

namespace Jetscape {

template<class T>
JetScapeReader<T>::JetScapeReader()
{
  VERBOSE(8);
  currentEvent=-1;
}

template<class T>
JetScapeReader<T>::~JetScapeReader()
{
  VERBOSE(8);
}

template<class T>
void JetScapeReader<T>::Clear()
{
  nodeVec.clear();edgeVec.clear();
  //pShower->clear();//pShower=nullptr; //check ...
  pShowers.clear();
}

template<class T>
void JetScapeReader<T>::AddNode(string s)
{
  string token; 
  //int counter=0;
  strT.set(s);

  vector<string> vS;
  
  while (!strT.done())
    {
      token = strT.next();
      vS.push_back(token);
      //++counter;
      //cout << token  << " ";
    }

  //cout<<endl;
  //cout<<vS.size()<<endl;
  nodeVec.push_back(pShower->new_vertex(make_shared<Vertex>(stod(vS[1]),stod(vS[2]),stod(vS[3]),stod(vS[4]))));
  //cout<<nodeVec.size()<<endl;
}

template<class T>
void JetScapeReader<T>::AddEdge(string s)
{
  if (nodeVec.size()>1)
    {
      string token; 
      //int counter=0;
      strT.set(s);
      
      vector<string> vS;
      
      while (!strT.done())
	{
	  token = strT.next();
	  vS.push_back(token);
	  //++counter;
	  //cout << token  << " ";
	}      
      //cout<<endl;
      
      pShower->new_parton(nodeVec[stoi(vS[0])],nodeVec[stoi(vS[1])],make_shared<Parton>(stoi(vS[2]),stoi(vS[3]),stoi(vS[4]),stod(vS[5]),stod(vS[6]),stod(vS[7]),stod(vS[8]))); // use different constructor wit true spatial posiiton ...
    }
  else
    WARN<<"Node vector not filled, can not add edges/partons!";
}

template<class T>
void JetScapeReader<T>::AddHadron(string s)
{
	string token;
        strT.set(s);
	token = strT.next();
	if(!token.compare("H"))
	{
		vector<string> vS;
		double x[4];
        	x[0]=x[1]=x[2]=x[3]=0.0;
		while (!strT.done())
        	{
			token = strT.next();
			vS.push_back(token);
		}
		hadrons.push_back(make_shared<Hadron>(stoi(vS[1]),stoi(vS[2]),stoi(vS[3]),stod(vS[4]),stod(vS[5]),stod(vS[6]),stod(vS[7]),x));
	}
}

template<class T>
void JetScapeReader<T>::Next()
{
  if (currentEvent>0)
    Clear();
  
  //ReadEvent(currentPos);
  //INFO<<"Read Events ...";
  string line;
  string token; 

  INFO<<"Current Event = "<<currentEvent;

  pShowers.push_back(make_shared<PartonShower>());
  pShower=pShowers[0];
  currentShower=1;
  
  int nodeZeroCounter=0;
  
  while (getline(inFile,line))
    {
      
      strT.set(line);
      
      if (!strT.isCommentEntry())
	{
	  if (strT.isEventEntry())
	    {
	      int newEvent=stoi(strT.next());		      
	      if (currentEvent!=newEvent && currentEvent>-1)
		{currentEvent++;break;}
	      
	      currentEvent=newEvent;	      	
	    }
	  else
	    {
	      if (strT.isGraphEntry())
		{
		  if (strT.isEdgeEntry())
		    {
		      //cout<<line<<endl;
		      AddEdge(line);
		    }
		  else		    
		    {
		      if (strT.isNodeZero())
			{
			  nodeZeroCounter++;
			  //cout<<nodeZeroCounter<<endl;
			  //cout<<line<<endl;
			  if (nodeZeroCounter>currentShower)
			    {
			      //pShower->clear();
			      //pShower=nullptr;
			      nodeVec.clear();edgeVec.clear();
			      pShowers.push_back(make_shared<PartonShower>());
			      pShower=pShowers.back();
			      currentShower++;
			    }
			}
		      AddNode(line);		     
		    }
		}
		else
		{
		 AddHadron(line);
		}
	    }
	}
    }
  //cout<<"There are "<<hadrons.size()<<" hadrons"<<endl; 
  if (Finished())
    currentEvent++;  
  
}

template<class T>
void JetScapeReader<T>::Init()
{
  VERBOSE(8)<<"Open Input File = "<<file_name_in;
  INFO<<"Open Input File = "<<file_name_in;
  
  inFile.open(file_name_in.c_str());
  
  if (!inFile.good())
    { WARN<<"Corrupt input file!"; exit(-1);}
  else
    INFO<<"File opened";

  currentEvent=0;
}

template class JetScapeReader<ifstream>;

template class JetScapeReader<igzstream>;

} // end namespace Jetscape
