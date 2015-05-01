//
//  Filename:     gauss_seidel.h
//  Programmer:   James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   7 - A Gauss-Seidel Iteration Class
//
//  Description:  This is the gauss-seidel iteration function class definition.
//

#ifndef __hw7__gauss_seidel__
#define __hw7__gauss_seidel__

#include "matrix.h"
#include "config.h"
#include "math.h"

template <class T>
class gauss_seidel
{
public:

	//Description: Perform Gaussian Elimination with Scaled Partial Pivoting
	//Pre:         matrix aA must be a (n x n) matrix
	//             vector aB must be of size n
	//Post:        vector aX will contain the size n solution to Ax = b
	bool operator()(vector<T>& aX, const matrix_base<T>& aA, const vector<T>& aB);
};

#include "gauss_seidel.hpp"

#endif
