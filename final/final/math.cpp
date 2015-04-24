//
//  Filename:     math.cpp
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  These are math convenience functions
//
//

#include "math.h"


double roundNumber(const double& aNumber, const double& aRoundTolerance)
{
	return floor(aNumber / aRoundTolerance + 0.5) * aRoundTolerance;
}

bool equivalent(const double& aLHS, const double& aRHS, const double& aEpsilon)
{
	return (fabs(aLHS - aRHS) < aEpsilon);
}

void seedRandom()
{
	srand(static_cast<unsigned int>(time(NULL)));
}