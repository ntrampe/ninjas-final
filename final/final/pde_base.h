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

	//// Known points of the pde

	//Description:  x lower bound function
	//Pre:          none
	//Post:         returns the known value of the pde at (bounds.x(), aY)
	virtual T xLower(T aY) const = 0;

	//Description:  x upper bound function
	//Pre:          none
	//Post:         returns the known value of the pde at (bounds.y(), aY)
	virtual T xUpper(T aY) const = 0;

	//Description:  y lower bound function
	//Pre:          none
	//Post:         returns the known value of the pde at (aX, bounds.x())
	virtual T yLower(T aX) const = 0;

	//Description:  y upper bound function
	//Pre:          none
	//Post:         returns the known value of the pde at (aX, bounds.y())
	virtual T yUpper(T aX) const = 0;

	// size of mesh
	size_t m_density;

	// bounds of the pde
	// x = lowerbound
	// y = upperbound
	point2d<T> m_bounds;

	// unknown points
	vector<point3d<T>> m_points;

public:

	//// Constructors

	//Description:  default constructor
	//Pre:          none
	//Post:         sets bounds to (0,1)
	pde_base();

	//Description:  bounds constructor
	//Pre:          none
	//Post:         sets bounds to aBounds
	pde_base(const size_t aN, point2d<T> aBounds);

	//Description:  bounds contructor
	//Pre:          none
	//Post:         sets bounds to (aLowerBound, aUpperBound)
	pde_base(const size_t aN, const T aLowerBound, const T aUpperBound);

	//Description:  destructor
	//Pre:          none
	//Post:         clears points
	virtual ~pde_base();


	//// Convienience


	//Description:  get density
	//Pre:          none
	//Post:         returns density
	size_t density() const;

	//Description:  get bounds
	//Pre:          none
	//Post:         returns bounds
	point2d<T> bounds() const;

	//Description:  get lower bound
	//Pre:          none
	//Post:         returns lower bound
	T lowerBound() const;

	//Description:  get upper bound
	//Pre:          none
	//Post:         returns upper bound
	T upperBound() const;
  
  
  //Description:  has the pde been solved
  //Pre:          none
  //Post:         returns true if any unknown points have been calculated
  bool solved() const;


	//Description:  set density
	//Pre:          none
	//Post:         sets density
	void setDensity(const size_t aN);
  
  //Description:  set bounds
  //Pre:          none
  //Post:         sets bounds
  void setBounds(const point2d<T> aBounds);

	//Description:  add unknown point
	//Pre:          none
	//Post:         adds (aX, aY, aZ) to points
	void addPoint(const T aX, const T aY, const T aZ);

	//Description:  add unknown point
	//Pre:          none
	//Post:         adds (aPoint.x(), aPoint.y(), aZ) to points
	void addPoint(const point2d<T>& aPoint, const T aZ);

	//Description:  clear points
	//Pre:          none
	//Post:         points cleared
	void clearPoints();

	//Description:  get matlab code to display graph of pde
	//Pre:          none
	//Post:         returns matlab code
	std::string matlabOutput(float aAnimationFactor = 0.0, const bool aDrawLines = false) const;

	//Description:  get unknown points
	//Pre:          none
	//Post:         returns unknown points
	std::string pointsOutput() const;

	//// Operators


	//Description:  projection operator
	//Pre:          (aX, aY) should be on the boundary
	//Post:         returns z value from boundary point
	T operator()(const T aX, const T aY) const;
};

#include "pde_base.hpp"

#endif /* defined(__final__pde_base__) */
