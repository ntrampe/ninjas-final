//
//  Filename:     cholesky.h
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   7 - A Parameterized Matrix Class and Cholesky Decomposition
//
//  Description:  This is the cholesky function class definition.
//                
//

#ifndef __final__cholesky__
#define __final__cholesky__

#include "matrix.h"
#include "config.h"
#include "math.h"

template <class T>
class cholesky
{
public:

	//Description: Perform Cholesky Decomposition
	//Pre:         matrix aA must be a (n x n) matrix
	//             vector aB must be of size n
	//			   division (/), *, +, and - operators must be defined for type T
	//Post:        vector aX will contain the size n solution to Ax = b
	bool operator()(vector<T>& aX, const matrix_base<T>& aA, const vector<T>& aB);
};

#include "cholesky.hpp"

#endif /* defined(__final__cholesky__) */
