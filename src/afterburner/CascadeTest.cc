#include "CascadeTest.h"
#include "JetScapeLogger.h"
#include "QueryHistory.h"

#include <thread>

using namespace Jetscape;

void CascadeTest::InitTask()
{
  JSINFO << "Cascade Test Init() dummy ...";
}

void CascadeTest::CalculateTime()
{
  VERBOSE(2) << "CascadeTest::CalculateTime() wait 100ms ... current time = " << GetModuleCurrentTime() << " Thread Id = " << std::this_thread::get_id();
  //create dummy hadrons ...
  for (int i=0; i<2; i++)
    hList.push_back(make_shared<Hadron>(0,211,0,i+1,0,0,i+1));

  std::this_thread::sleep_for(std::chrono::milliseconds(100));
}

void CascadeTest::ExecTime()
{
  vector<any> eLossHistories = QueryHistory::Instance()->GetHistoryFromModules("JetEnergyLoss");

  if (GetMainClock()->GetCurrentTime() < 2) {

    cout << "CascadeTest::ExecTime(): Current Main Clock Time = " << GetMainClock()->GetCurrentTime() << endl;
    cout << "CascadeTest::ExecTime(): Print Histories via vector<any> eLossHistories = QueryHistory::Instance()->GetHistoryFromModules(\"JetEnergyLoss\")" << endl;

    for (auto mHist : eLossHistories)
    {
      any_cast<std::shared_ptr<PartonShower>>(mHist)->PrintEdges(false);
    }
  }
}

void CascadeTest::InitPerEvent()
{

}

void CascadeTest::FinishPerEvent()
{
  hList.clear();
}

any CascadeTest::GetHistory()
{
  return any(hList);
}