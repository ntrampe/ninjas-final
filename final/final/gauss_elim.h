//
//  Filename:     gauss_elim.h
//  Programmer:   James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  This is the gaussian elimination function class definition.
//                This function class performs Gaussian Elimination
//                with Scaled Partial Pivoting
//

#ifndef __hw7__gauss_elim__
#define __hw7__gauss_elim__

#include "matrix.h"
#include "config.h"
#include "math.h"

template <class T>
class gauss_elim
{
public:

	//Description: Perform Gaussian Elimination with Scaled Partial Pivoting
	//Pre:         matrix aA must be a (n x n) matrix
	//             vector aB must be of size n
	//Post:        vector aX will contain the size n solution to Ax = b
	bool operator()(vector<T>& aX, const matrix_base<T>& aA, const vector<T>& aB);
};

#include "gauss_elim.hpp"

#endif
