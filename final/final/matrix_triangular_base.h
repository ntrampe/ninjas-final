//
//  Filename:     matrix_triangular_base.h
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  This is the matrix_triangular_base abstract base class
//                The purpose of this class is to provide an easy way
//                to create a child matrix class without redefining
//                a large amount of implementation.
//

#ifndef __hw6__matrix_triangular_base__
#define __hw6__matrix_triangular_base__

#include "matrix_base.h"

template <class T>
class matrix_triangular_base : public matrix_base<T>
{
protected:

	// size of matrix
	size_t m_size;

	// the element returned for an element located outside
	// of the triangular area
	T m_outsideElement;

	//Description: Convert (aRow, aColumn) to the corresponding index in m_data
	//Pre:         aRow and aColumns should be within the dimensions of the matrix
	//Post:        aIndex is set the the index containing the (aRow, aColumn) element
	virtual void convertCoordinatesToIndex(size_t& aIndex, const size_t aRow, const size_t aColumn) const = 0;

	//Description: Convert aIndex to the corresponding (aRow, aColumn)
	//Pre:         aIndex should be within the bounds of the data array (area of matrix)
	//Post:        aRow and aColumn are set to the (aRow, aColumn) containing the
	//             m_data[aIndex] element
	virtual void convertIndexToCoordinates(size_t& aRow, size_t& aColumn, const size_t aIndex) const = 0;

	//Description: Determine if an element is in the triangular area
	//Pre:         aRow and aColumn should be within the dimensions of the matrix
	//Post:        returns true if an element is in the triangular area
	virtual bool withinTriangle(const size_t aRow, const size_t aColumn) const = 0;

	//Description: Retrieve an element to modify
	//Pre:         aRow and aColumn must be within the matrix bounds
	//Post:        returns a reference to the element at (aRow, aColumn)
	//             if (aRow, aColumn) is outside the triangular area,
	//             m_outsideElement is returned
	virtual T& at(const size_t aRow, const size_t aColumn);

	//Description: Retrieve an element to read only
	//Pre:         aRow and aColumn must be within the matrix bounds
	//Post:        returns a copy of the element at (aRow, aColumn)
	//             if (aRow, aColumn) is outside the triangular area,
	//             m_outsideElement is returned
	virtual T at(const size_t aRow, const size_t aColumn) const;

	//Description: Allocate and initialize memory for data array
	//Pre:         if aMemorySize == 0, then memorySize() will be used
	//Post:        creates a new m_data array and sets the elements to T()
	virtual void initMatrix(const size_t aMemorySize = 0);

	//Description: Set up member variables in terms of rows and columns
	//Pre:         none
	//Post:        sets size to min of aRow and aColumn
	virtual void setupMatrix(const size_t aRow, const size_t aColumn);

public:


	//// Constructors / Destructor


	//Description: Default constructor
	//Pre:         none
	//Post:        creates a (1 x 1) empty matrix
	matrix_triangular_base();

	//Description: Dimension contructor
	//Pre:         aRows and aColumns should be greater than 0
	//Post:        creates a (aRows x aColumns) matrix
	matrix_triangular_base(const size_t aSize);

	//Description: Copy constructor
	//Pre:         none
	//Post:        copies aCopy's dimensions and data over to matrix
	matrix_triangular_base(const matrix_triangular_base<T>& aCopy);

	//Description: Destructor
	//Pre:         none
	//Post:        deallocates data array
	~matrix_triangular_base();


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
	virtual std::string name() const = 0;


	//// Column / Row Manipulation


	//Description: Replace (aRow, aVector[i]) with the contents of aVector
	//Pre:         aVector's size must be equal to the matrix's columns
	//             aRow must be within matrix dimensions
	//Post:        the matrix will have the contents of aVector at aRow
	//             only traverses triangular elements
	virtual void replaceVectorAtRow(const vector<T>& aVector, const size_t aRow) = 0;

	//Description: Replace (aVector[i], aColumn) with the contents of aVector
	//Pre:         aVector's size must be equal to the matrix's rows
	//             aColumn must be within matrix dimensions
	//Post:        the matrix will have the contents of aVector at aColumn
	//             only traverses triangular elements
	virtual void replaceVectorAtColumn(const vector<T>& aVector, const size_t aColumn) = 0;


	//Description: Solve the matrix as a system of equations
	//Pre:         vector aB's size must be equal to the matrix rows
	//Post:        aX will contain the solution of Ax = b, x = aX
	virtual bool solveMatrix(const vector<T>& aB, vector<T>& aX) = 0;

	//Description: Assignment operator
	//Pre:         class used in template needs to overload = operator
	//Post:        sets matrix dimensions and data to aRHS'
	matrix_triangular_base<T>& operator=(const matrix_triangular_base<T>& aRHS);
};

#include "matrix_triangular_base.hpp"

#endif /* defined(__matrix_triangular_base__) */
