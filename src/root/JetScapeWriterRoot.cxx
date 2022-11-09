// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// jetscape writer root class

#include "JetScapeWriterRoot.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include "MakeUniqueHelper.h"

#include "JetScapeRootParton.h"
#include "JetScapeRootVertex.h"
#include "JetScapeRootPartonShower.h"

//#include "TXMLEngine.h"

using namespace std;

//Just for DEBBUG here ...
/*
void DisplayNode(TXMLEngine &xml, XMLNodePointer_t node, Int_t level)
{
   // this function display all accessible information about xml node and its children
   printf("%*c node: %s\n", level, ' ', xml.GetNodeName(node));
   // display namespace
   XMLNsPointer_t ns = xml.GetNS(node);
   if (ns != 0)
      printf("%*c namespace: %s refer: %s\n", level + 2, ' ', xml.GetNSName(ns), xml.GetNSReference(ns));
   // display attributes
   XMLAttrPointer_t attr = xml.GetFirstAttr(node);
   while (attr != 0) {
      printf("%*c attr: %s value: %s\n", level + 2, ' ', xml.GetAttrName(attr), xml.GetAttrValue(attr));
      attr = xml.GetNextAttr(attr);
   }
   // display content (if exists)
   const char *content = xml.GetNodeContent(node);
   if (content != 0)
      printf("%*c cont: %s\n", level + 2, ' ', content);
   // display all child nodes
   XMLNodePointer_t child = xml.GetChild(node);
   while (child != 0) {
      DisplayNode(xml, child, level + 2);
      child = xml.GetNext(child);
   }
}
*/

ClassImp(JetScapeWriterRoot)

JetScapeWriterRoot::JetScapeWriterRoot(string m_file_name_out)
{
  SetOutputFileName(m_file_name_out);
}

JetScapeWriterRoot::~JetScapeWriterRoot()
{
  VERBOSE(8);
  if (GetActive())
      Close();
}

void JetScapeWriterRoot::WriteEvent()
{
  JSDEBUG<< GetCurrentEvent() << " Event";
  //Write(to_string(GetCurrentEvent()) + " Event");
  //event=make_unique<JetScapeRootEvent>();
  event->SetEventNumber(GetCurrentEvent());

}

/*
void JetScapeWriterRoot::Write(weak_ptr<Parton> p)
{
  //output_file<<*p.lock()<<endl;
}

void JetScapeWriterRoot::Write(weak_ptr<Vertex> v)
{
  //output_file<<*v.lock()<<endl;
}
*/

void JetScapeWriterRoot::Write(weak_ptr<PartonShower> pS)
{

  JSDEBUG<<"Write PartonShower in ROOT ...";

  PartonShower::node_iterator nIt,nEnd;

  auto mS=make_unique<JetScapeRootPartonShower>();

  for (nIt = pS.lock()->nodes_begin(), nEnd = pS.lock()->nodes_end(); nIt != nEnd; ++nIt)
    {
      auto mv=make_unique<JetScapeRootVertex>(*pS.lock()->GetVertex(*nIt));
      mv->SetNodeId(nIt->id());

      mS->AddVertex(move(mv));

      //DEBUG:
      //mv->Write("mv");
    }

  PartonShower::edge_iterator eIt,eEnd;

  //DEBUG:
  //cout<<pS.lock()-> GetNumberOfPartons() <<endl;

  for (eIt = pS.lock()->edges_begin(), eEnd = pS.lock()->edges_end(); eIt != eEnd; ++eIt)
    {
      auto mp=make_unique<JetScapeRootParton>(*pS.lock()->GetParton(*eIt));
      mp->SetSourceNode(eIt->source().id());mp->SetTargetNode(eIt->target().id());

      mS->AddParton(move(mp));
      //DEBUG:
      //mp->Write("mp");
    }

  //mS->Write("ms");
  event->AddPartonShower(move(mS));
}

void JetScapeWriterRoot::Init()
{
   if (GetActive())
     {
       JSINFO<<"JetScape Root Writer initialized with output file = "<<GetOutputFileName();
       output_file=make_unique<TFile>(GetOutputFileName().c_str(),"RECREATE");
       //output_file->SetCompressionLevel(9);

       event=new JetScapeRootEvent();
       evTree=make_shared<TTree>("eventTree","JetScape Event Tree");
       evTree->Branch("Event",&event);

       //Write Init Informations, like XML and ... to file ...
       //WriteInitFileXML();
     }
}

void JetScapeWriterRoot::Exec()
{
  JSINFO<<"Run JetScapeWriterRoot: Write event # "<<GetCurrentEvent()<<" ...";

  if (GetActive())
    WriteEvent();
}

void JetScapeWriterRoot::Clear()
{
  JSDEBUG<<"Fill ROOT Tree ...";
  evTree->Fill();
  delete event;

  event=new JetScapeRootEvent();
}

void JetScapeWriterRoot::WriteInitFileXML()
{
  /*
  DEBUG<<"Write XML to output file. XML file = "<<JetScapeXML::Instance()->GetXMLFileName();
  tinyxml2::XMLPrinter printer;
  JetScapeXML::Instance()->GetXMLDocument().Print(&printer);
  //cout<<printer.CStr()<<endl;
  WriteComment("Init XML file used : "+JetScapeXML::Instance()->GetXMLFileName());
  output_file<<printer.CStr();
  */

  JSINFO<<"Write XML to output file. XML file = "<<JetScapeXML::Instance()->GetXMLUserFileName();
  tinyxml2::XMLPrinter printer;
  JetScapeXML::Instance()->GetXMLDocumentUser().Print(&printer);
  //cout<<printer.CStr()<<endl;

  //Remark: Not ideal ...
  TNamed mXMLfile("jetscape_init",printer.CStr());

  //cout<<mXMLfile<<endl;
  mXMLfile.Write();

  /*
  TXMLEngine xml;
  // Now try to parse xml file
  // Only file with restricted xml syntax are supported
  XMLDocPointer_t xmldoc = xml.ParseFile((JetScapeXML::Instance()->GetXMLFileName()).c_str());
  if (!xmldoc) return;
  // take access to main node
  XMLNodePointer_t mainnode = xml.DocGetRootElement(xmldoc);
  // display recursively all nodes and subnodes
  //DisplayNode(xml, mainnode, 1);
  // Release memory before exit

  xml.Write("jetscape_init");
  xml.FreeDoc(xmldoc);
  */
}
