//
//  Filename:     matrix_symmetrical.h
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  This is the matrix_symmetrical class definition.
//                The matrix_symmetrical template class represents an (n x n) symmetrical matrix
//

#ifndef __hw6__matrix_symmetrical__
#define __hw6__matrix_symmetrical__

#include "matrix_base.h"

template <class T>
class matrix_symmetrical : public matrix_base<T>
{
protected:

	// size of matrix
	size_t m_size;

	//Description: Convert (aRow, aColumn) to the corresponding index in m_data
	//Pre:         aRow and aColumns should be within the dimensions of the matrix
	//Post:        aIndex is set the the index containing the (aRow, aColumn) element
	virtual void convertCoordinatesToIndex(size_t& aIndex, const size_t aRow, const size_t aColumn) const;

	//Description: Convert aIndex to the corresponding (aRow, aColumn)
	//Pre:         aIndex should be within the bounds of the data array (area of matrix)
	//Post:        aRow and aColumn are set to the (aRow, aColumn) containing the
	//             m_data[aIndex] element
	virtual void convertIndexToCoordinates(size_t& aRow, size_t& aColumn, const size_t aIndex) const;

	//Description: Retrieve an element to modify
	//Pre:         aRow and aColumn must be within the matrix bounds
	//Post:        returns a reference to the element at (aRow, aColumn)
	//             if (aRow, aColumn) is outside the triangular area,
	//             the symmetrical is returned
	virtual T& at(const size_t aRow, const size_t aColumn);

	//Description: Retrieve an element to read only
	//Pre:         aRow and aColumn must be within the matrix bounds
	//Post:        returns a copy of the element at (aRow, aColumn)
	//             if (aRow, aColumn) is outside the triangular area,
	//             the symmetrical is returned
	virtual T at(const size_t aRow, const size_t aColumn) const;

	//Description: Set up member variables in terms of rows and columns
	//Pre:         none
	//Post:        sets size to min of aRow and aColumn
	virtual void setupMatrix(const size_t aRow, const size_t aColumn);

public:


	//// Constructors / Destructor


	//Description: Default constructor
	//Pre:         none
	//Post:        creates a (1 x 1) empty matrix
	matrix_symmetrical();

	//Description: Dimension contructor
	//Pre:         aRows and aColumns should be greater than 0
	//Post:        creates a (aRows x aColumns) matrix
	matrix_symmetrical(const size_t aSize);

	//Description: Copy constructor
	//Pre:         none
	//Post:        copies aCopy's dimensions and data over to matrix
	matrix_symmetrical(const matrix_symmetrical<T>& aCopy);

	//Description: Destructor
	//Pre:         none
	//Post:        deallocates data array
	virtual ~matrix_symmetrical();


	//// Convenience


	//Description: Retrieve rows
	//Pre:         none
	//Post:        returns number of rows in matrix
	virtual size_t rows() const;

	//Description: Retrieve columns
	//Pre:         none
	//Post:        returns number of columns in matrix
	virtual size_t columns() const;

	//Description: Retrieve memorySize (((1 + size) * size) / 2)
	//Pre:         none
	//Post:        returns the size of the data array
	virtual size_t memorySize() const;

	//Description: Retrieve size
	//Pre:         none
	//Post:        returns number of columns/rows in matrix
	size_t size() const;

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
	//             uses Cholesky Decomposition
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
	matrix_symmetrical<T>& operator=(const matrix_base<T>& aRHS);

	//Description: Negation operator
	//Pre:         class used in template needs to overload - operator
	//Post:        returns a matrix that is the negative of aRHS
	template <class U>
	friend matrix_symmetrical<U> operator-(const matrix_symmetrical<U>& aRHS);

	//Description: Default addition operator for any derived matrix class
	//             including matrix_symmetrical
	//Pre:         class used in template needs to overload + operator
	//Post:        returns a matrix that is the addition of aLHS and aRHS
	template <class U>
	friend matrix_symmetrical<U> operator+(const matrix_symmetrical<U>& aLHS, const matrix_symmetrical<U>& aRHS);

	//Description: Subtraction operator
	//Pre:         class used in template needs to overload - operator
	//Post:        returns a matrix that is the subtraction of aLHS and aRHS
	template <class U>
	friend matrix_symmetrical<U> operator-(const matrix_symmetrical<U>& aLHS, const matrix_symmetrical<U>& aRHS);

	//Description: Scalar Multiplication operator
	//Pre:         class used in template needs to overload * (double) operator
	//Post:        returns a matrix where all elements of aRHS are multiplied by aLHS
	template <class U>
	friend matrix_symmetrical<U> operator*(const double& aLHS, const matrix_symmetrical<U>& aRHS);

	//Description: Scalar Multiplication operator
	//Pre:         class used in template needs to overload * (double) operator
	//Post:        returns a matrix where all elements of aLHS are multiplied by aRHS
	template <class U>
	friend matrix_symmetrical<U> operator*(const matrix_symmetrical<U>& aLHS, const double& aRHS);

	//Description: Vector Multiplication operator
	//Pre:         class used in template needs to overload * (vector) operator
	//Post:        returns a matrix that is the multiplication of aLHS and aRHS
	template <class U>
	friend vector<U> operator*(const vector<U>& aLHS, const matrix_symmetrical<U>& aRHS);

	//Description: Vector Multiplication operator
	//Pre:         class used in template needs to overload * (vector) operator
	//Post:        returns a matrix that is the multiplication of aRHS and aLHS
	template <class U>
	friend vector<U> operator*(const matrix_symmetrical<U>& aLHS, const vector<U>& aRHS);
};

#include "matrix_symmetrical.hpp"

#endif /* defined(__matrix_symmetrical__) */
