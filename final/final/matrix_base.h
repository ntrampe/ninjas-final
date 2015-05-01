//
//  Filename:     matrix_base.h
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  This is the matrix_base abstract base class
//                The purpose of this class is to provide an easy way
//                to create a child matrix class without redefining
//                a large amount of implementation.
//
//                The idea is that every matrix class will have a
//                data array:
//
//                  m_data = T[0 .. n]    where n = memorySize()
//
//                and two one-to-one functions to convert between the
//                data array index and element location:
//
//                  (i, j)  ->  [n]
//                  [n]     ->  (i, j)
//                                        where n >= 0 && n < memorySize()
//                                              i >= 0 && i < this->rows()
//                                              j >= 0 && j < this->columns()
//
//                This way, we can move a significant amount of functionality
//                into this base class, rather than implementing it in every
//                child class as redundant code. And, we still keep memory
//                efficiency:
//
//                  ex)
//
//                    n x m arbitrary matrix:
//                      sizeof(matrix<T>(n, m)) = sizeof(T)*(n * m)
//
//                    n x n upper triangular matrix:
//                      sizeof(matrix_triangular_upper<T>(n)) = sizeof(T)*(((1 + n) * n) / 2 + 1)
//

#ifndef __hw6__matrix_base__
#define __hw6__matrix_base__

#include <iostream>
#include <iomanip>
#include <string>

#include "math.h"
#include "vector.h"
#include "point2d.h"

template <class T>
class matrix_base
{
protected:

	// data array
	T* m_data = nullptr;

	//Description: Convert (aRow, aColumn) to the corresponding index in m_data
	//Pre:         aRow and aColumns should be within the dimensions of the matrix
	//Post:        aIndex is set the the index containing the (aRow, aColumn) element
	virtual void convertCoordinatesToIndex(size_t& aIndex, const size_t aRow, const size_t aColumn) const = 0;

	//Description: Convert aIndex to the corresponding (aRow, aColumn)
	//Pre:         aIndex should be within the bounds of the data array
	//Post:        aRow and aColumn are set to the (aRow, aColumn) containing the
	//             m_data[aIndex] element
	virtual void convertIndexToCoordinates(size_t& aRow, size_t& aColumn, const size_t aIndex) const = 0;

	//Description: Determine if aRow and aColumns are within the matrix dimensions
	//Pre:         none
	//Post:        returns true if aRow < rows() and aColumn < aColumns
	virtual bool withinDimensions(const size_t aRow, const size_t aColumn) const;

	//Description: Retrieve an element to modify
	//Pre:         aRow and aColumn must be within the matrix bounds
	//Post:        returns a reference to the element at (aRow, aColumn)
	virtual T& at(const size_t aRow, const size_t aColumn);

	//Description: Retrieve an element to read only
	//Pre:         aRow and aColumn must be within the matrix bounds
	//Post:        returns a copy of the element at (aRow, aColumn)
	virtual T at(const size_t aRow, const size_t aColumn) const;

	//Description: Set up subclass member variables in terms of rows and columns
	//Pre:         see derived class
	//Post:        see derived class
	virtual void setupMatrix(const size_t aRow, const size_t aColumn) = 0;

	//Description: Allocate and initialize memory for data array
	//Pre:         if aMemorySize == 0, then memorySize() will be used
	//Post:        creates a new m_data array and sets the elements to T()
	virtual void initMatrix(const size_t aMemorySize = 0);

	//Description: Init matrix and copy over respective elements
	//Pre:         memorySize() should be set up correctly before
	//             e.g.) set rows and columns before calling
	//Post:        creates a new m_data array and sets the elements to aMatrix's
	virtual void copyMatrix(const matrix_base<T>& aMatrix);

	//Description: Retrieve longest element size in matrix (for output)
	//Pre:         none
	//Post:        returns longest element size in matrix
	size_t longestElement() const;

public:

	//// Destructor

	//Description: Destructor
	//Pre:         none
	//Post:        deallocates data array
	virtual ~matrix_base();

	//// Convenience

	//Description: Retrieve rows
	//Pre:         none
	//Post:        returns number of rows in matrix
	virtual size_t rows() const = 0;

	//Description: Retrieve columns
	//Pre:         none
	//Post:        returns number of columns in matrix
	virtual size_t columns() const = 0;

	//Description: Retrieve memory size
	//Pre:         none
	//Post:        returns the length of the data array for
	//             a particular subclass
	virtual size_t memorySize() const = 0;

	//Description: Name of matrix
	//Pre:         none
	//Post:        returns user-defined name of matrix class
	virtual std::string name() const = 0;

	//Description: Description of matrix
	//Pre:         none
	//Post:        returns user-defined description of matrix class
	virtual std::string description() const;
  
  size_t lengthOfDiagonal(const size_t aRow) const;

	//Description: Clear the matrix
	//Pre:         none
	//Post:        sets all data to the default value
	virtual void clear(const T aDefault = T());

	//Description: Randomize the values in the matrix (for testing purposes, mainly)
	//Pre:         none
	//Post:        sets all the elements to random values no larger than aMaximum
	virtual void randomize(const unsigned int aMaximum = 10);

	//Description: Order the values in the matrix (for testing purposes, mainly)
	//Pre:         none
	//Post:        sets all the elements to incremental values starting at aStart
	virtual void order(const T aStart = T());

  //Description: Prints matrix memory in order
  //Pre:         none
  //Post:        outputs all the elements in index order
	void printMemory() const;

	//Description: If any element is near zero, set it to exactly zero
	//Pre:         templated class should
	//Post:        sets all the elements to close to zero equal to zero
	virtual void roundToZero();

	//Description: Resize the matrix while preserving the data
	//Pre:         none
	//Post:        sets rows to aRows and columns to aColumns
	//             each element stays in its old (row, column) location
	void resize(const size_t aRows, const size_t aColumns);


	//Description: Determine equality by rows, columns and elements,
	//             but not data size / derived class
	//             this allows different derived classes to be equal to each other
	//
	//             Always (rows * columns) time, use == operator for efficiency
	//
	//Pre:         none
	//Post:        returns true if all elements are equal
	bool isEqualTo(const matrix_base<T>& aMatrix) const;


	//// Column / Row Manipulation


	//Description: Retrieve a vector at a row
	//Pre:         aRow must be within matrix dimensions
	//Post:        returns a vector of the elements at aRow
	vector<T> vectorAtRow(const size_t aRow) const;

	//Description: Retrieve a vector at a column
	//Pre:         aColumn must be within matrix dimensions
	//Post:        returns a vector of the elements at aColumn
	vector<T> vectorAtColumn(const size_t aColumn) const;

	//Description: Replace (aRow, aVector[i]) with the contents of aVector
	//Pre:         aVector's size must be equal to the matrix's columns
	//             aRow must be within matrix dimensions
	//Post:        the matrix will have the contents of aVector at aRow
	virtual void replaceVectorAtRow(const vector<T>& aVector, const size_t aRow);

	//Description: Replace (aVector[i], aColumn) with the contents of aVector
	//Pre:         aVector's size must be equal to the matrix's rows
	//             aColumn must be within matrix dimensions
	//Post:        the matrix will have the contents of aVector at aColumn
	virtual void replaceVectorAtColumn(const vector<T>& aVector, const size_t aColumn);


	//// Row operations


	//Description: Solve the matrix as a system of equations
	//Pre:         vector aB's size must be equal to the matrix rows
	//Post:        aX will contain the solution of Ax = b, x = aX
	virtual bool solveMatrix(const vector<T>& aB, vector<T>& aX) = 0;

	//Description: Determine if the matrix is in row echelon form
	//Pre:         none
	//Post:        returns true if the pivots are zero and it's an upper
	//             triagular matrix
	bool isRowEchelonForm() const;

	//Description: Switch two rows
	//Pre:         aSource and aDestination must be within the matrix bounds
	//Post:        contents of the row at aSource are switched with the contents
	//             of the row at aDestination
	void switchRows(const size_t aSource, const size_t aDestination);

	//Description: Scale a row
	//Pre:         aSource and aDestination must be within the matrix bounds
	//Post:        contents of the row at aDestination contain the contents of
	//             aSource * aScalar
	void scaleRow(const size_t aSource, const double aScalar, const size_t aDestination);

	//Description: Add two rows
	//Pre:         aLHS, aRHS and aDestination must be within the matrix bounds
	//Post:        contents of the row at aDestination contain the contents of
	//             the row at aLHS + the contents of the row at aRHS
	void addRow(const size_t aLHS, const size_t aRHS, const size_t aDestination);


	//// Indexing Operators


	//Description: Retrieve an element to modify
	//Pre:         aRow and aColumn must be within the matrix bounds
	//Post:        returns a reference to the element at (aRow, aColumn)
	T& operator()(const size_t aRow, const size_t aColumn);

	//Description: Retrieve an element to read only
	//Pre:         aRow and aColumn must be within the matrix bounds
	//Post:        returns a copy of the element at (aRow, aColumn)
	T operator()(const size_t aRow, const size_t aColumn) const;

	//Description: Retrieve an element to modify
	//Pre:         aPoint must be within the matrix bounds
	//Post:        returns a reference to the element at (aPoint.x(), aPoint.y())
	T& operator()(const point2d<size_t>& aPoint);

	//Description: Retrieve an element to read only
	//Pre:         aPoint must be within the matrix bounds
	//Post:        returns a copy of the element at (aPoint.x(), aPoint.y())
	T operator()(const point2d<size_t>& aPoint) const;

	//// Operators

	//Description: Assignment operator
	//Pre:         class used in template needs to overload = operator
	//Post:        sets matrix dimensions and data to aRHS'
	//
	//NOTE:        Each derived class will only need to forward their assignment operator's
	//             functionality to this base class
	//             Unique member variable setup will go into setupMatrix
	//             This forces a creator of a derived class to set up their member variables
	//             before calling copy matrix... unless they don't forward the functionality...
	matrix_base<T>& operator=(const matrix_base<T>& aRHS);

	//Description: Comparison Operator
	//Pre:         class used in template needs to overload != operator
	//Post:        returns whether or not all elements in aLHS are equal to
	//             all elements in aRHS
	template <class U>
	friend bool operator==(const matrix_base<U>& aLHS, const matrix_base<U>& aRHS);

	//Description: Not comparison operator
	//Pre:         none
	//Post:        returns whether or not any elements in aLHS are not equal to
	//             any elements in aRHS
	template <class U>
	friend bool operator!=(const matrix_base<U>& aLHS, const matrix_base<U>& aRHS);

	//Description: Output Stream operator
	//Pre:         class used in template needs to overload << operator
	//Post:        outputs the matrix to the console
	template <class U>
	friend std::ostream& operator<<(std::ostream& aOutput, const matrix_base<U>& aMatrix);

	//Description: Input Stream operator
	//Pre:         stream must be in the proper format
	//              0 .. m
	//              . .  .
	//              .  . .
	//              n .. n,m
	//             where n = the number of rows, m = the number of columns
	//             class used in template needs to overload >> operator
	//             aMatrix should be initialized with the rows and columns
	//             desired for input
	//Post:        aMatrix is set to the contents aInput
	template <class U>
	friend std::istream& operator>>(std::istream& aInput, matrix_base<U>& aMatrix);
};

#include "matrix_base.hpp"

#endif /* defined(__hw6__matrix_base__) */
