//
//  Filename:     runtime.cpp
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  
//

#include "runtime.h"


void runtime::begin()
{
  gettimeofday(&m_begin, NULL);
}


void runtime::end()
{
  gettimeofday(&m_end, NULL);
}


double runtime::elapsed()
{
  //convert to microseconds
  double elapsedTime = (m_end.tv_sec - m_begin.tv_sec) * 1000.0;
  elapsedTime += (m_end.tv_usec - m_begin.tv_usec) / 1000.0;
  return elapsedTime;
}
