//
//  Filename:     runtime.h
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  Class used to get measure the runtime of different solving methods
//

#ifndef __Sorting__nt_runtime__
#define __Sorting__nt_runtime__

#include <iostream>
#include <sys/time.h> // for gettimeofday()

class runtime
{
private:
  timeval m_begin;
  timeval m_end;
  
public:
  //Pre: none
  //Post: sets the beginning time to the current time
  //Description: Set the starting time
  void begin();
  
  //Pre: none
  //Post: sets the ending time to the current time
  //Description: Set the ending time
  void end();
  
  //Pre: the beginning and ending time should already be set
  //Post: returns the time difference between begin and end
  //Description: Get the elapsed time
  double elapsed();
};

#endif /* defined(__Sorting__nt_runtime__) */
