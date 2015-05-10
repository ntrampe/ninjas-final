//
//  Filename:     matrix_triangular_lower.h
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  This class represents a lower triangular
//                matrix.
//

#ifndef __final__matrix_triangular_lower__
#define __final__matrix_triangular_lower__

#include "matrix_triangular_base.h"

template <class T>
class matrix_triangular_lower : public matrix_triangular_base<T>
{
private:

	//Description: Convert (aRow, aColumn) to the corresponding index in m_data
	//Pre:         aRow and aColumns should be within the dimensions of the matrix
	//Post:        aIndex is set the the index containing the (aRow, aColumn) element
	virtual void convertCoordinatesToIndex(size_t& aIndex, const size_t aRow, const size_t aColumn) const;

	//Description: Convert aIndex to the corresponding (aRow, aColumn)
	//Pre:         aIndex should be within the bounds of the data array
	//Post:        aRow and aColumn are set to the (aRow, aColumn) containing the
	//             m_data[aIndex] element
	virtual void convertIndexToCoordinates(size_t& aRow, size_t& aColumn, const size_t aIndex) const;

	//Description: Determine if an element is in the triangular area
	//Pre:         aRow and aColumn should be within the dimensions of the matrix
	//Post:        returns true if an element is in the triangular area
	virtual bool withinTriangle(const size_t aRow, const size_t aColumn) const;

public:

	//// Constructors / Destructor


	//Description: Default constructor
	//Pre:         none
	//Post:        creates a (1 x 1) empty matrix
	matrix_triangular_lower() : matrix_triangular_base<T>() {}

	//Description: Dimension contructor
	//Pre:         aRows and aColumns should be greater than 0
	//Post:        creates a (aRows x aColumns) matrix filled with aDefault
	matrix_triangular_lower(const size_t aSize) : matrix_triangular_base<T>(aSize) {}

	//Description: Copy constructor
	//Pre:         none
	//Post:        copies aCopy's dimensions and data over to matrix
	matrix_triangular_lower(const matrix_triangular_lower<T>& aCopy) : matrix_triangular_base<T>(aCopy) {}


	//// Convienience


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


	//// Operators


	//Description: Base assignment operator
	//Pre:         class used in template needs to overload = operator
	//Post:        sets matrix dimensions and data to aRHS'
	matrix_triangular_lower<T>& operator=(const matrix_base<T>& aRHS);

	//Description: Negation operator
	//Pre:         class used in template needs to overload - operator
	//Post:        returns a matrix that is the negative of aRHS
	template <class U>
	friend matrix_triangular_lower<U> operator-(const matrix_triangular_lower<U>& aRHS);

	//Description: Addition operator
	//Pre:         class used in template needs to overload + operator
	//Post:        returns a matrix that is the addition of aLHS and aRHS
	template <class U>
	friend matrix_triangular_lower<U> operator+(const matrix_triangular_lower<U>& aLHS, const matrix_triangular_lower<U>& aRHS);

	//Description: Subtraction operator
	//Pre:         class used in template needs to overload - operator
	//Post:        returns a matrix that is the subtraction of aLHS and aRHS
	template <class U>
	friend matrix_triangular_lower<U> operator-(const matrix_triangular_lower<U>& aLHS, const matrix_triangular_lower<U>& aRHS);

	//Description: Scalar Multiplication operator
	//Pre:         class used in template needs to overload * (double) operator
	//Post:        returns a matrix where all elements of aRHS are multiplied by aLHS
	template <class U>
	friend matrix_triangular_lower<U> operator*(const double& aLHS, const matrix_triangular_lower<U>& aRHS);

	//Description: Scalar Multiplication operator
	//Pre:         class used in template needs to overload * (double) operator
	//Post:        returns a matrix where all elements of aLHS are multiplied by aRHS
	template <class U>
	friend matrix_triangular_lower<U> operator*(const matrix_triangular_lower<U>& aLHS, const double& aRHS);

	//Description: Matrix Multiplication operator
	//Pre:         class used in template needs to overload * and + operator
	//Post:        returns a matrix that is the addition of aLHS and aRHS
	template <class U>
	friend matrix_triangular_lower<U> operator*(const matrix_triangular_lower<U>& aLHS, const matrix_triangular_lower<U>& aRHS);

	//Description: Vector Multiplication operator
	//Pre:         class used in template needs to overload * (vector) operator
	//Post:        returns a matrix that is the multiplication of aLHS and aRHS
	template <class U>
	friend vector<U> operator*(const vector<U>& aLHS, const matrix_triangular_lower<U>& aRHS);

	//Description: Vector Multiplication operator
	//Pre:         class used in template needs to overload * (vector) operator
	//Post:        returns a matrix that is the multiplication of aRHS and aLHS
	template <class U>
	friend vector<U> operator*(const matrix_triangular_lower<U>& aLHS, const vector<U>& aRHS);
};

#include "matrix_triangular_lower.hpp"

#endif /* defined(__final__matrix_triangular_lower__) */
