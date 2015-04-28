//
//  Filename:     gaussSPP.h
//  Programmer:   Nicholas Trampe
//  Class:        CS 5201 - Clayton Price
//  Assignment:   4 - A Parameterized Matrix Class and Gaussian Elimination
//
//  Description:  This is the gaussSPP function class definition.
//                This function class performs Gaussian Elimination
//                with Scaled Partial Pivoting
//

#ifndef __hw4__gaussSPP__
#define __hw4__gaussSPP__

#include "matrix.h"
#include "config.h"
#include "math.h"

template <class T>
class gaussSPP
{
public:

	//Description: Perform Gaussian Elimination with Scaled Partial Pivoting
	//Pre:         matrix aA must be a (n x n) matrix
	//             vector aB must be of size n
	//Post:        vector aX will contain the size n solution to Ax = b
	bool operator()(vector<T>& aX, const matrix<T>& aA, const vector<T>& aB);
};

#include "nt_gaussSPP.hpp"

#endif /* defined(__hw4__gaussSPP__) */
