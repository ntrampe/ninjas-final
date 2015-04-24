//
//  Filename:     matrix_diagonal.h
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  This is the matrix_diagonal class definition.
//                The matrix_diagonal template class represents an (n x n) diagonal matrix
//

#ifndef __hw6__matrix_diagonal__
#define __hw6__matrix_diagonal__

#include "matrix_banded.h"

template <class T>
class matrix_diagonal : public matrix_banded<T>
{
protected:

	//Description: Convert (aRow, aColumn) to the corresponding index in m_data
	//Pre:         aRow and aColumns should be within the dimensions of the matrix
	//Post:        aIndex is set the the index containing the (aRow, aColumn) element
	virtual void convertCoordinatesToIndex(size_t& aIndex, const size_t aRow, const size_t aColumn) const;

	//Description: Convert aIndex to the corresponding (aRow, aColumn)
	//Pre:         aIndex should be within the bounds of the data array (area of matrix)
	//Post:        aRow and aColumn are set to the (aRow, aColumn) containing the
	//             m_data[aIndex] element
	virtual void convertIndexToCoordinates(size_t& aRow, size_t& aColumn, const size_t aIndex) const;

public:


	//// Constructors / Destructor


	//Description: Default constructor
	//Pre:         none
	//Post:        creates a (1 x 1) empty matrix
	matrix_diagonal() : matrix_banded<T>(1, 0) {}

	//Description: Dimension contructor
	//Pre:         aRows and aColumns should be greater than 0
	//Post:        creates a (aRows x aColumns) matrix filled with aDefault
	matrix_diagonal(const size_t aSize) : matrix_banded<T>(aSize, 0) {}

	//Description: Copy constructor
	//Pre:         none
	//Post:        copies aCopy's dimensions and data over to matrix
	matrix_diagonal(const matrix_diagonal<T>& aCopy);

	//Description: Base copy constructor
	//Pre:         none
	//Post:        copies aCopy's dimensions and data over to matrix
	matrix_diagonal(const matrix_base<T>& aMatrix);


	//// Convenience

	//Description: Name of matrix
	//Pre:         none
	//Post:        returns user-defined name of matrix class
	virtual std::string name() const;


	//// Column / Row Manipulation

	//Description: Replace (aRow, aVector[i]) with the contents of aVector
	//Pre:         aVector's size must be equal to the matrix's columns
	//             aRow must be within matrix dimensions
	//Post:        the matrix will have the contents of aVector at aRow
	//             only traverses triangular elements
	virtual void replaceVectorAtRow(const vector<T>& aVector, const size_t aRow);

	//Description: Replace (aVector[i], aColumn) with the contents of aVector
	//Pre:         aVector's size must be equal to the matrix's rows
	//             aColumn must be within matrix dimensions
	//Post:        the matrix will have the contents of aVector at aColumn
	//             only traverses triangular elements
	virtual void replaceVectorAtColumn(const vector<T>& aVector, const size_t aColumn);


	//Description: Solve the matrix as a system of equations
	//             scales pivots to one
	//Pre:         vector aB's size must be equal to the matrix rows
	//Post:        aX will contain the solution of Ax = b, x = aX
	virtual bool solveMatrix(const vector<T>& aB, vector<T>& aX);


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
	matrix_diagonal<T>& operator=(const matrix_base<T>& aRHS);

	//Description: Matrix Multiplication operator
	//Pre:         class used in template needs to overload * and + operator
	//Post:        returns a matrix that is the addition of aLHS and aRHS
	template <class U>
	friend matrix_diagonal<U> operator*(const matrix_diagonal<U>& aLHS, const matrix_diagonal<U>& aRHS);

	//Description: Vector Multiplication operator
	//Pre:         class used in template needs to overload * (vector) operator
	//Post:        returns a matrix that is the multiplication of aLHS and aRHS
	template <class U>
	friend vector<U> operator*(const vector<U>& aLHS, const matrix_diagonal<U>& aRHS);

	//Description: Vector Multiplication operator
	//Pre:         class used in template needs to overload * (vector) operator
	//Post:        returns a matrix that is the multiplication of aRHS and aLHS
	template <class U>
	friend vector<U> operator*(const matrix_diagonal<U>& aLHS, const vector<U>& aRHS);
};

#include "matrix_diagonal.hpp"

#endif /* defined(__hw6__matrix_diagonal__) */
