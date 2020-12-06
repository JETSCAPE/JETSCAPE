#include "BrickTest.h"
#include <thread>

using namespace Jetscape;

void BrickTest::CalculateTime()
{
    VERBOSE(2)<<"BrickTest::CalculateTime() wait 100ms ... current time = "<<GetModuleCurrentTime()<<" Thread Id = "<<std::this_thread::get_id();
    EvolveHydro();
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
}

void BrickTest::ExecTime()
{

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