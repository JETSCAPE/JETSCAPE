#include "BrickTest.h"
#include "QueryHistory.h"
#include <thread>

using namespace Jetscape;

void BrickTest::CalculateTime()
{
	VERBOSE(2) << "BrickTest::CalculateTime() wait 100ms ... current time = " << GetModuleCurrentTime() << " Thread Id = " << std::this_thread::get_id();
	EvolveHydro();
	std::this_thread::sleep_for(std::chrono::milliseconds(100));
}

void BrickTest::ExecTime()
{
	if (GetMainClock()->GetCurrentTime() < 2) {

    cout << "BrickTest::ExecTime(): Current Main Clock Time = " << GetMainClock()->GetCurrentTime() << endl;
    cout << "BrickTest::ExecTime(): Print Histories via vector<any> eLossHistories = QueryHistory::Instance()->GetHistoryFromModules(\"CascadeTest\")" << endl;

    cout<<"BrickTest::ExecTime(): # of Hadrons in cascade = "<<any_cast<std::vector<std::shared_ptr<Hadron>>>(QueryHistory::Instance()->GetHistoryFromModule("CascadeTest")).size()<<endl;
  }
}

void BrickTest::InitPerEvent()
{

}

void BrickTest::FinishPerEvent()
{

}

any BrickTest::GetHistory()
{
	return 0;
}