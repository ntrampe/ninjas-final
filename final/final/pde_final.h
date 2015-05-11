//
//  Filename:     pde_final.h
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  This is the pde class definition.
//                The pde class represents a partial differential equation
//

#ifndef __final__pde_final__
#define __final__pde_final__

#include "pde_base.h"

template <class T>
class pde_final : public pde_base<T>
{
private:

	// boundary functions
	virtual T xLower(T aY) const;
	virtual T xUpper(T aY) const;
	virtual T yLower(T aX) const;
	virtual T yUpper(T aX) const;

public:

	//Description:  default constructor
	//Pre:          none
	//Post:         sets bounds to (0,1)
	pde_final() : pde_base<T>() {}

	//Description:  bounds constructor
	//Pre:          none
	//Post:         sets bounds to aBounds
	pde_final(const size_t aN, point2d<T> aBounds) : pde_base<T>(aN, aBounds) {}

	//Description:  bounds contructor
	//Pre:          none
	//Post:         sets bounds to (aLowerBound, aUpperBound)
	pde_final(const size_t aN, const T aLowerBound, const T aUpperBound) : pde_base<T>(aN, aLowerBound, aUpperBound) {}
};

#include "pde_final.hpp"

#endif /* defined(__final__pde_final__) */
