//
//  Filename:     point.h
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  This is the point class definition.
//                The point template class represents a (x, y) point
//                in two-dimensional space
//                Or for this project, a location in a matrix
//

#ifndef __hw6__poi_
#define __hw6__poi_

template <class T>
class point
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
	//Post:        sets point location to default
	point() : m_x(T()), m_y(T()) {}

	//Description: Point constructor
	//Pre:         none
	//Post:        sets point location to (aX, aY)
	point(T aX, T aY) : m_x(aX), m_y(aY) {}

	//Description: Destructor
	//Pre:         none
	//Post:        deconstructs object
	~point();

	//// Convenience

	//Description: Set point location
	//Pre:         none
	//Post:        sets point location to (aX, aY)
	void set(const T aX, const T aY);

	//Description: Set x location
	//Pre:         none
	//Post:        sets point location to (aX, m_y)
	void setX(const T aX);

	//Description: Set y location
	//Pre:         none
	//Post:        sets point location to (m_y, aY)
	void setY(const T aY);

	//Description: Get x
	//Pre:         none
	//Post:        returns x location
	T x() const;

	//Description: Get y
	//Pre:         none
	//Post:        returns y location
	T y() const;
  
  //Pre:         none
  //             class used in template needs to overload << operator
  //Post:        outputs aVector elements to aOutput
  //Description: output vector data
  template <class U>
  friend std::ostream& operator<<(std::ostream& aOutput, const point<U>& aPoint);
};

#include "point.hpp"

#endif /* defined(__hw6__poi_) */
