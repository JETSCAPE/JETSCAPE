#ifndef Stopwatch_header
#define Stopwatch_header

#include <ctime>

class Stopwatch
{
  private:
    time_t start, end;
  public:
    Stopwatch() {start=clock(); end=0;}
    void tic() {start=clock();}
    void toc() {end=clock();}
    double takeTime() {return ((double)(end - start)) / CLOCKS_PER_SEC;}
};

#endif

/*-----------------------------------------------------------------------
  Usage:
  Declare a class as:
    Stopwatch sw;
  Then the time a piece of code takes can be recorded as:
    sw.tic();
    **** code ****
    sw.toc();
  And the result can be outputted using:
    cout << sw.takeTime() << endl;
-----------------------------------------------------------------------*/
