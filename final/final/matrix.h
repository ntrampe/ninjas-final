//
//  Filename:     matrix.h
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  This is the matrix class definition.
//                The matrix template class represents an (n x m) matrix,
//                where:
//
//                A =
//
//                  A(0,0) .. A(0,m)
//
//                  .   .     .
//
//                  .      .  .
//
//                  A(n,0) .. A(n,m)
//
//                n = the number of rows
//                m = the number of columns
//
//                A(n, m) = the element at the nth row and the mth column
//

#ifndef __hw6__matrix__
#define __hw6__matrix__

#include "matrix_base.h"

template <class T>
class matrix : public matrix_base<T>
{
protected:

	// number of rows
	size_t m_rows;

	// number of columns
	size_t m_columns;

	//Description: Convert (aRow, aColumn) to the corresponding index in m_data
	//Pre:         aRow and aColumns should be within the dimensions of the matrix
	//Post:        aIndex is set the the index containing the (aRow, aColumn) element
	virtual void convertCoordinatesToIndex(size_t& aIndex, const size_t aRow, const size_t aColumn) const;

	//Description: Convert aIndex to the corresponding (aRow, aColumn)
	//Pre:         aIndex should be within the bounds of the data array (area of matrix)
	//Post:        aRow and aColumn are set to the (aRow, aColumn) containing the
	//             m_data[aIndex] element
	virtual void convertIndexToCoordinates(size_t& aRow, size_t& aColumn, const size_t aIndex) const;

	//Description: Set up member variables in terms of rows and columns
	//Pre:         none
	//Post:        sets rows and columns to aRow and aColumn
	virtual void setupMatrix(const size_t aRow, const size_t aColumn);

public:


	//// Constructors / Destructor


	//Description: Default constructor
	//Pre:         none
	//Post:        creates a (1 x 1) empty matrix
	matrix();

	//Description: Dimension contructor
	//Pre:         aRows and aColumns should be greater than 0
	//Post:        creates a (aRows x aColumns) matrix filled with aDefault
	matrix(const size_t aRows, const size_t aColumns);

	//Description: Copy constructor
	//Pre:         none
	//Post:        copies aCopy's dimensions and data over to matrix
	matrix(const matrix<T>& aCopy);

	//Description: Base copy constructor
	//Pre:         none
	//Post:        copies aCopy's dimensions and data over to matrix
	matrix(const matrix_base<T>& aMatrix);

	//Description: Vector constructor
	//Pre:         none
	//Post:        create a (aVector.size() x 1) column matrix filled with
	//             the elements of aVector
	matrix(const vector<T>& aVector);

	//Description: Destructor
	//Pre:         none
	//Post:        deallocates data array and sets dimensions to 0
	virtual ~matrix();


	//// Convenience


	//Description: Retrieve rows
	//Pre:         none
	//Post:        returns number of rows in matrix
	virtual size_t rows() const;

	//Description: Retrieve columns
	//Pre:         none
	//Post:        returns number of columns in matrix
	virtual size_t columns() const;

	//Description: Retrieve memorySize (row * columns)
	//Pre:         none
	//Post:        returns the size of the data array
	virtual size_t memorySize() const;

	//Description: Name of matrix
	//Pre:         none
	//Post:        returns user-defined name of matrix class
	virtual std::string name() const;


	//// Column / Row Manipulation

	//Description: Insert aVector after aRow
	//Pre:         aVector's size must be equal to the matrix's columns
	//Post:        the matrix will have at least one more column
	//             the matrix will have the contents of aVector at aRow
	void insertVectorAtRow(const vector<T>& aVector, const size_t aRow);

	//Description: Insert aVector after aColumn
	//Pre:         aVector's size must be equal to the matrix's rows
	//Post:        the matrix will have at least one more row
	//             the matrix will have the contents of aVector at aColumn
	void insertVectorAtColumn(const vector<T>& aVector, const size_t aColumn);

	//Description: Get a matrix within the bounds specified
	//Pre:         matrix must be a square matrix
	//             start and end rows/columns must be within the matrix bounds
	//Post:        returns a matrix containing the elements between the specified bounds
	matrix<T> getSubMatrix(const size_t aStartRow, const size_t aStartColumn, const size_t aEndRow, const size_t aEndColumn) const;

	//Description: Get the matrix's determinant
	//Pre:         the matrix must be square
	//Post:        returns the matrix determinant
	double determinant() const;


	//Description: Solve the matrix as a system of equations
	//             uses Gaussian Elimination with scaled partial pivoting
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
	matrix<T>& operator=(const matrix_base<T>& aRHS);

	//Description: Negation operator
	//Pre:         class used in template needs to overload - operator
	//Post:        returns a matrix that is the negative of aRHS
	template <class U>
	friend matrix<U> operator-(const matrix_base<U>& aRHS);

	//Description: Transpose operator
	//Pre:         none
	//Post:        returns the transpose of aRHS
	template <class U>
	friend matrix<U> operator~(const matrix_base<U>& aRHS);

	//Description: Default addition operator for any derived matrix class
	//             including matrix
	//             Other, more efficient operators are defined in the child classes
	//             themselves
	//Pre:         class used in template needs to overload + operator
	//Post:        returns a matrix that is the addition of aLHS and aRHS
	template <class U>
	friend matrix<U> operator+(const matrix_base<U>& aLHS, const matrix_base<U>& aRHS);

	//Description: Subtraction operator
	//Pre:         class used in template needs to overload - operator
	//Post:        returns a matrix that is the subtraction of aLHS and aRHS
	template <class U>
	friend matrix<U> operator-(const matrix_base<U>& aLHS, const matrix_base<U>& aRHS);

	//Description: Matrix Multiplication operator
	//Pre:         class used in template needs to overload * and + operator
	//Post:        returns a matrix that is the addition of aLHS and aRHS
	template <class U>
	friend matrix<U> operator*(const matrix_base<U>& aLHS, const matrix_base<U>& aRHS);

	//Description: Scalar Multiplication operator
	//Pre:         class used in template needs to overload * (double) operator
	//Post:        returns a matrix where all elements of aRHS are multiplied by aLHS
	template <class U>
	friend matrix<U> operator*(const double& aLHS, const matrix_base<U>& aRHS);

	//Description: Scalar Multiplication operator
	//Pre:         class used in template needs to overload * (double) operator
	//Post:        returns a matrix where all elements of aLHS are multiplied by aRHS
	template <class U>
	friend matrix<U> operator*(const matrix_base<U>& aLHS, const double& aRHS);

	//Description: Vector Multiplication operator
	//Pre:         class used in template needs to overload * (vector) operator
	//Post:        returns a matrix that is the multiplication of aLHS and aRHS
	template <class U>
	friend vector<U> operator*(const vector<U>& aLHS, const matrix_base<U>& aRHS);

	//Description: Vector Multiplication operator
	//Pre:         class used in template needs to overload * (vector) operator
	//Post:        returns a matrix that is the multiplication of aRHS and aLHS
	template <class U>
	friend vector<U> operator*(const matrix_base<U>& aLHS, const vector<U>& aRHS);
};

#include "matrix.hpp"

#endif /* defined(__hw6__matrix__) */
