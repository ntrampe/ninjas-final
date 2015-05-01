//
//  Filename:     vector.h
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  This is the vector class definition.
//                The vector template class represents a dynamic array
//

#ifndef __hw6__vector__
#define __hw6__vector__

#include <stdio.h>
#include <cstdlib>
#include <stdexcept>
#include <iostream>
#include <cstring>
#include <sstream>
#include <cmath>
#include <limits>

template <class T>
class vector
{
private:

	//// Variables

	// size of elements in vector
	size_t m_size;

	// size of allocated memory in terms of number of elements
	size_t m_maxSize;

	// data pointer
	T* m_data;

	//// Memory Management Functions

	//Pre:         none
	//Post:        allocates a larger array and copies data if size
	//             is 0 or one less than max size
	//Description: grow memory array
	void addMemory();

	//Pre:         none
	//Post:        allocates a smaller array and copies data if size
	//             is the next smallest power of 2 from max size
	//Description: shrink memory array
	void removeMemory();

	//Pre:         none
	//Post:        sets max size to the closest, but larger power of two
	//             to aMemory, allocates a new array of that size and
	//             copies data over, if applicable
	//Description: allocate new data array
	void setMemory(const size_t aMemory);

public:

	//// Constructors / Destructor

	//Pre:         none
	//Post:        initializes size and data to 0/NULL
	//Description: default constructor
	vector() : m_size(0), m_maxSize(0), m_data(NULL) {}

	//Pre:         none
	//Post:        initializes size and data to copy's
	//             size and data
	//Description: copy constructor
	vector(const vector<T>& aCopy);

	//Pre:         none
	//Post:        deallocates data
	//Description: destructor
	~vector();

	//// Core Functions

	//Pre:         none
	//Post:        sets size and data to 0/NULL
	//Description: clear vector elements
	void clear();

	//Pre:         none
	//Post:        sets data/memory size to aQuantity
	//             if aUpdateSize is true, the size of the
	//             vector will be modified
	//             otherwise, the only the memory will be allocated
	//Description: reserve memory for a vector
	void reserve(const size_t aQuantity, const bool aUpdateSize = false);

	//Pre:         none
	//Post:        inserts aElement into the last index
	//Description: append element to the back of the vector
	void push_back(const T& aElement);

	//Pre:         none
	//Post:        removes aElement from the vector, if it exists
	//             shifts remaining elements back
	//Description: remove an element from a vector
	void remove(const T& aElement);

	//Pre:         vector must not be empty
	//Post:        remove element from the back of the vector
	//Description: remove back element
	void pop_back();

	//Pre:         none
	//Post:        returns the size of the elements in the vector
	//Description: size of vector elements
	size_t size() const;

	//Pre:         none
	//Post:        returns the size of all allocated memory
	//Description: size of allocated memory (in element quantity)
	size_t maxSize() const;

	//// Convenience

	//Pre:         class used in template needs to be able to convert to a double
	//Post:        returns whether or not a vector is a multiple
	//Description: size of allocated memory (in element quantity)
	bool isMultiple(const vector<T>& aRHS) const;

	//// Operators

	//Pre:         aIndex must be >= 0 and < m_size
	//Post:        returns a reference to the aIndexth element
	//             user can modify the element directly
	//Description: retrieve and modify the element at aIndex
	T& operator[](size_t aIndex);

	//Pre:         aIndex must be >= 0 and < m_size
	//Post:        returns a constant reference to the aIndexth element
	//Description: retrieve the element at aIndex
	const T& operator[](size_t aIndex) const;

	//Pre:         none
	//             class used in template needs to overload << operator
	//Post:        outputs aVector elements to aOutput
	//Description: output vector data
	template <class U>
	friend std::ostream& operator<<(std::ostream& aOutput, const vector<U>& aVector);

	//Pre:         aInput must be in the form 0 .. n
	//             aVector should be initialized with the size
	//             desired for input
	//Post:        populates aVector with contents of aInput
	//Description: input vector data
	template <class U>
	friend std::istream& operator>>(std::istream& aInput, vector<U>& aVector);

	//Pre:         class used in template needs to overload + operator
	//Post:        returns the addition of aLHS and aRHS by adding repective elements together
	//Description: add two vectors
	template <class U>
	friend vector<U> operator+(const vector<U>& aLHS, const vector<U>& aRHS);

	//Pre:         class used in template needs to overload - operator
	//Post:        returns the subtraction of aLHS and aRHS by subtracting repective elements together
	//Description: subtract two vectors
	template <class U>
	friend vector<U> operator-(const vector<U>& aLHS, const vector<U>& aRHS);

	//Pre:         class used in template needs to overload * operator
	//Post:        returns the dot product of aLHS and aRHS
	//Description: dot product of two vectors
	template <class U>
	friend U operator*(const vector<U>& aLHS, const vector<U>& aRHS);

	//Pre:         class used in template needs to overload * with a double on the left
	//Post:        returns a vector containing the multiplication of every element
	//             in aRHS with aLHS
	//Description: multiply a vector by a constant
	template <class U>
	friend vector<U> operator*(const double aLHS, const vector<U>& aRHS);

	//Pre:         class used in template needs to overload * with a double on the right
	//Post:        returns a vector containing the multiplication of every element
	//             in aLHS with aRHS
	//Description: multiply a vector by a constant
	template <class U>
	friend vector<U> operator*(const vector<U>& aLHS, const double aRHS);

	//Pre:         class used in template needs to overload negation operator
	//Post:        returns the negation of aVector
	//Description: negate vector
	template <class U>
	friend vector<U> operator-(const vector<U>& aVector);

	//Pre:         class used in template needs to overload != operator
	//Post:        returns whether or not all elements in aLHS are equal to
	//             all elements in aRHS
	//Description:
	template <class U>
	friend bool operator==(const vector<U>& aLHS, const vector<U>& aRHS);

	//Pre:         none
	//Post:        returns whether or not any elements in aLHS are not equal to
	//             any elements in aRHS
	//Description:
	template <class U>
	friend bool operator!=(const vector<U>& aLHS, const vector<U>& aRHS);

	//Pre:         class used in template needs to overload = operator
	//Post:        sets vector size and data to aRHS's size and data
	//Description: assignment operator
	vector<T>& operator=(const vector<T>& aRHS);
};

#include "vector.hpp"

#endif /* defined(__hw6__vector__) */
