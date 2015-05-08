//
//  Filename:     pde_base.h
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  This is the pde class definition.
//                The pde class represents a partial differential equation
//

#ifndef __final__pde_base__
#define __final__pde_base__

#include <stdio.h>
#include "vector.h"
#include "point2d.h"
#include "point3d.h"
#include "math.h"

template <class T>
class pde_base
{
private:
  
  virtual T xLower(T aY) const = 0;
  virtual T xUpper(T aY) const = 0;
  virtual T yLower(T aX) const = 0;
  virtual T yUpper(T aX) const = 0;
  
  point2d<T> m_bounds;
  vector<point3d<T>> m_points;
  
public:
  
  pde_base();
  pde_base(point2d<T> aBounds);
  pde_base(const T aLowerBound, const T aUpperBound);
  virtual ~pde_base();
  
  point2d<T> bounds() const;
  T lowerBound() const;
  T upperBound() const;
  
  void addPoint(const T aX, const T aY, const T aZ);
  void addPoint(const point2d<T>& aPoint, const T aZ);
  
  void addKnownPoint(const T aX, const T aY);
  void addKnownPoint(const point2d<T>& aPoint);
  
  void clearPoints();
  
  std::string matlabOutput(float aAnimationFactor = 0.0, const bool aDrawLines = false) const;
  
  T operator()(const T aX, const T aY) const;
};

#include "pde_base.hpp"

#endif /* defined(__final__pde_base__) */
