//
//  Filename:     pde_notes.h
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  This is the pde class definition.
//                The pde class represents a partial differential equation
//

#ifndef __final__pde_notes__
#define __final__pde_notes__

#include "pde_base.h"

template <class T>
class pde_notes : public pde_base<T>
{
protected:
  
  virtual T xLower(T aY) const;
  virtual T xUpper(T aY) const;
  virtual T yLower(T aX) const;
  virtual T yUpper(T aX) const;
  
public:
  
  pde_notes() : pde_base<T>() {}
  pde_notes(point2d<T> aBounds) : pde_base<T>(aBounds) {}
  pde_notes(const T aLowerBound, const T aUpperBound) : pde_base<T>(aLowerBound, aUpperBound) {}
};

#include "pde_notes.hpp"

#endif /* defined(__final__pde_notes__) */
