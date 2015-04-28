//
//  Filename:     pdeBounds.h
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  This is the pde class definition.
//                The pde class represents a partial differential equation
//

#ifndef __final__pde__
#define __final__pde__

#include <stdio.h>
#include "point.h"

template <class T>
class pdeBounds
{
private:
  T (*m_xLower)(T);
  T (*m_xUpper)(T);
  T (*m_yLower)(T);
  T (*m_yUpper)(T);
  
  point<T> m_bounds;
  
public:
  
  pdeBounds();
  pdeBounds(T (*aXLower)(T), T (*aXUpper)(T), T (*aYLower)(T), T (*aYUpper)(T), point<T> aBounds);
  ~pdeBounds();
  
  T operator()(const T aX, const T aY);
};

#include "pdeBounds.hpp"

#endif /* defined(__final__pde__) */
