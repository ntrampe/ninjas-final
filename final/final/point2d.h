//
//  Filename:     point2d.h
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  This is the point2d class definition.
//                The point2d template class represents a (x, y) point2d
//                in two-dimensional space
//                Or for this project, a location in a matrix
//

#ifndef __hw6__poi_
#define __hw6__poi_

#include <iomanip>

template <class T>
class point2d
{
private:

	// x location
	T m_x;

	// y location
	T m_y;

public:

	//// Constructors / Destructor

	//Description: Default constructor
	//Pre:         none
	//Post:        sets point2d location to default
	point2d() : m_x(T()), m_y(T()) {}

	//Description: point2d constructor
	//Pre:         none
	//Post:        sets point2d location to (aX, aY)
	point2d(T aX, T aY) : m_x(aX), m_y(aY) {}

	//Description: Destructor
	//Pre:         none
	//Post:        deconstructs object
	~point2d();

	//// Convenience

	//Description: Set point2d location
	//Pre:         none
	//Post:        sets point2d location to (aX, aY)
	void set(const T aX, const T aY);

	//Description: Set x location
	//Pre:         none
	//Post:        sets point2d location to (aX, m_y)
	void setX(const T aX);

	//Description: Set y location
	//Pre:         none
	//Post:        sets point2d location to (m_y, aY)
	void setY(const T aY);

	//Description: Get x
	//Pre:         none
	//Post:        returns x location
	T x() const;

	//Description: Get y
	//Pre:         none
	//Post:        returns y location
	T y() const;
  
  point2d<T>& operator=(const point2d<T>& aRHS);
  template <class U>
  friend bool operator==(const point2d<U>& aLHS, const point2d<U>& aRHS);
  template <class U>
  friend bool operator!=(const point2d<U>& aLHS, const point2d<U>& aRHS);
  
  //Pre:         none
  //             class used in template needs to overload << operator
  //Post:        outputs aVector elements to aOutput
  //Description: output vector data
  template <class U>
  friend std::ostream& operator<<(std::ostream& aOutput, const point2d<U>& aPoint);
};

#include "point2d.hpp"

#endif /* defined(__hw6__poi_) */
