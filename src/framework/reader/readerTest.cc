// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// Reader test (focus on graph)

#include <iostream>
#include <fstream>
#include <memory>
#include <chrono>
#include <thread>

#include "gzstream.h"
#include "PartonShower.h"
#include "JetScapeLogger.h"
#include "JetScapeReader.h"
#include "JetScapeBanner.h"

#include <GTL/dfs.h>

using namespace std;

// -------------------------------------

// Forward declaration
void Show();
void AnalyzeGraph(shared_ptr<PartonShower> mS);

// -------------------------------------

// Create a pdf of the shower graph:
// Use with graphviz (on Mac: brew install graphviz --with-app)
// in shell: dot GVfile.gv -Tpdf -o outputPDF.pdf
// [or you can also use the GraphViz app for Mac Os X (in the "cellar" of homebrew)]

// -------------------------------------

int main(int argc, char** argv)
{
  JetScapeLogger::Instance()->SetDebug(true);
  JetScapeLogger::Instance()->SetRemark(true);
  //SetVerboseLevel (9 a lot of additional debug output ...)
  //If you want to suppress it: use SetVerboseLevle(0) 
  JetScapeLogger::Instance()->SetVerboseLevel(9);
  
  cout<<endl;
  Show();

  vector<shared_ptr<PartonShower>> mShowers;

  //Directly with template: provide the relevant stream
  //auto reader=make_shared<JetScapeReader<ifstream>>("test_out.dat");
  //auto reader=make_shared<JetScapeReader<igzstream>>("test_out.dat.gz");
  
  // Hide Template (see class declarations in reader/JetScapeReader.h) ...
  auto reader=make_shared<JetScapeReaderAscii>("test_out.dat");
  //auto reader=make_shared<JetScapeReaderAsciiGZ>("test_out.dat.gz");

  // reads in multiple events and multiple shower per event
  // commentend out so that you get the dot graph file for the first shower in the first event
  // (add in and the file gets overriden)
  //while (!reader->Finished())
    {
      reader->Next();
      
      cout<<"Analyze current event = "<<reader->GetCurrentEvent()<<endl;
      mShowers=reader->GetPartonShowers();     
      
      for (int i=0;i<mShowers.size();i++)
	{
	  cout<<" Analyze parton shower = "<<i<<endl;
	  mShowers[i]->PrintNodes();
	  mShowers[i]->PrintEdges();	  
	  AnalyzeGraph(mShowers[i]);

	  if (i==0)	 
	    mShowers[i]->SaveAsGV("my_test.gv");	    

	  // wait for 5s
	  //std::this_thread::sleep_for(std::chrono::milliseconds(5000));  
	}
    }
 
}

// -------------------------------------

void AnalyzeGraph(shared_ptr<PartonShower> mS)
{
  INFO<<"Some GTL graph/shower analysis/dfs search output:";

  // quick and dirty ...
  dfs search;
  search.calc_comp_num(true);
  search.scan_whole_graph(true);
  search.start_node();// defaulted to first node ...
  search.run(*mS);

  cout<<endl;
  cout<<"DFS graph search feature from GTL:"<<endl;
  cout<<"Number of Nodes reached from node 0 = "<< search.number_of_reached_nodes () <<endl;
  cout<<"Node/Vertex ordering result from DFS:"<<endl;
  dfs::dfs_iterator itt2, endt2;
  for (itt2 = search.begin(), endt2=search.end(); itt2 !=endt2; ++itt2)
    {
      cout<<*itt2<<" ";//<<"="<<search.dfs_num(*itt2)<<" ";
    }
  cout<<endl;
  cout<<"Edge/Parton ordering result from DFS:"<<endl;
  dfs::tree_edges_iterator itt, endt;
  for (itt = search.tree_edges_begin(), endt=search.tree_edges_end(); itt !=endt; ++itt)
    {
        cout<<*itt;//<<endl;
    }
  cout<<endl;
  
  dfs::roots_iterator itt3, endt3;
  cout<<"List of root nodes found in graph/shower : ";
  for (itt3 = search.roots_begin(), endt3=search.roots_end(); itt3 !=endt3; ++itt3)
    {
      cout<<**itt3;
    }
  cout<<endl;
  cout<<endl;
}


// -------------------------------------

void Show()
{
  show_jetscape_banner();
  INFO_NICE;
  INFO_NICE<<"------------------------------------";
  INFO_NICE<<"| Reader Test JetScape Framework ... |";
  INFO_NICE<<"------------------------------------";
  INFO_NICE;
}
