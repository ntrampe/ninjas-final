//
//  Filename:     gauss_seidel.h
//  Programmer:   James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  This is the gauss-seidel iteration function class definition.
//

#ifndef __final__gauss_seidel__
#define __final__gauss_seidel__

#include "matrix.h"
#include "config.h"
#include "math.h"

template <class T>
class gauss_seidel
{
public:

	//Description: Default constructor for gauss_seidel class
	//Pre:         None
	//Post:        Creates an instance of gauss_seidel with error tolerance of 0.0000001
	gauss_seidel();

	//Description: Constructor for gauss_seidel class
	//Pre:         None
	//Post:        Creates an instance of gauss_seidel with error tolerance equalt to parameter
	gauss_seidel(const double aErrorTolerance);

	//Description: Destructor for gauss_seidel
	//Pre:         gauss_seidel has not been destroyed already
	//Post:        gauss_seidel is destroyed and memory is freed
	~gauss_seidel();

	//Description: Perform Gauss-Seidel Iteration
	//Pre:         matrix aA must be a (n x n) matrix
	//             vector aB must be of size n
	//			   division (/), *, +=, and - operators must be defined for type T
	//Post:        vector aX will contain the size n solution to Ax = b
	bool operator()(vector<T>& aX, const matrix_base<T>& aA, const vector<T>& aB);

private:
	double m_error_tol;

};

#include "gauss_seidel.hpp"

#endif
