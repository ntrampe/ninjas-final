//
//  Filename:     matrix_tridiagonal.h
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  This is the matrix_tridiagonal class definition.
//                The matrix_tridiagonal template class represents an (n x n) tridiagonal matrix
//

#ifndef __final__matrix_tridiagonal__
#define __final__matrix_tridiagonal__

#include "matrix_banded.h"

template <class T>
class matrix_tridiagonal : public matrix_banded<T>
{
protected:


public:


	//// Constructors / Destructor


	//Description: Default constructor
	//Pre:         none
	//Post:        creates a (1 x 1) empty matrix
	matrix_tridiagonal() : matrix_banded<T>(3, 1) {}

	//Description: Dimension contructor
	//Pre:         aSize should be greater than 0
	//Post:        creates a (aRows x aColumns) matrix filled with aDefault
	matrix_tridiagonal(const size_t aSize) : matrix_banded<T>(aSize, 1) {}

	//Description: Copy constructor
	//Pre:         none
	//Post:        copies aCopy's dimensions and data over to matrix
	matrix_tridiagonal(const matrix_tridiagonal<T>& aCopy);

	//Description: Base copy constructor
	//Pre:         none
	//Post:        copies aCopy's dimensions and data over to matrix
	matrix_tridiagonal(const matrix_base<T>& aMatrix);


	//// Convenience

	//Description: Name of matrix
	//Pre:         none
	//Post:        returns user-defined name of matrix class
	virtual std::string name() const;


	//// Operators

	// NOTE:  I use the base class in all the operations because
	//        the operations are the same. They would not be less efficient
	//        for other derived classes because those derived classes
	//        implement other class specific operators.
	//        These are for either dense matrices or a combination
	//        of any other derived class (e.g dense = upper + lower).


	//Description: Base assignment operator
	//Pre:         class used in template needs to overload = operator
	//Post:        sets matrix dimensions and data to aRHS'
	matrix_tridiagonal<T>& operator=(const matrix_base<T>& aRHS);
};

#include "matrix_tridiagonal.hpp"

#endif /* defined(__final__matrix_tridiagonal__) */
