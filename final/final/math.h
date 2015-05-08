//
//  Filename:     math.h
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  These are math convenience functions
//
//

#ifndef __final__math__
#define __final__math__

#include <iostream>
#include <cmath>
#include "config.h"

// what place to round numbers in roundNumber
const double DEFAULT_ROUND_TOLERANCE = 0.000001;

// epsilon in equivalent
const double DEFAULT_EQUIVALENCE_EPSILON = 0.00001;

//Pre:         none
//Post:        returns aNumber rounded to the aRoundTolerance place
//Description: Round a number
double roundNumber(const double& aNumber, const double& aRoundTolerance = DEFAULT_ROUND_TOLERANCE);

//Pre:         none
//Post:        returns whether or not aLHS is aEpsilon away from aRHS
//Description: Determine if two numbers are close to each other
bool equivalent(const double& aLHS, const double& aRHS, const double& aEpsilon = DEFAULT_EQUIVALENCE_EPSILON);

//Pre:         none
//Post:        seeds random with time at null
//Description: Seed random
void seedRandom();

#endif /* defined(__final__math__) */
