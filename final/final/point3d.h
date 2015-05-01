//
//  Filename:     point3d.h
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  This is the point3d class definition.
//                The point3d template class represents a (x, y) point3d
//                in two-dimensional space
//                Or for this project, a location in a matrix
//

#ifndef __hw6__point3d__
#define __hw6__point3d__

#include "point2d.h"

template <class T>
class point3d : public point2d<T>
{
private:

	// z location
	T m_z;

public:

	//// Constructors / Destructor

	//Description: Default constructor
	//Pre:         none
	//Post:        sets point3d location to default
	point3d() : point2d<T>(T(), T()), m_z(T()) {}

	//Description: point3d constructor
	//Pre:         none
	//Post:        sets point3d location to (aX, aY, aZ)
	point3d(T aX, T aY, T aZ) : point2d<T>(aX, aY), m_z(aZ) {}

	//Description: Destructor
	//Pre:         none
	//Post:        deconstructs object
	~point3d();

	//// Convenience

	//Description: Set point3d location
	//Pre:         none
	//Post:        sets point3d location to (aX, aY, aZ)
	void set(const T aX, const T aY, const T aZ);
  
  //Description: Set z location
  //Pre:         none
  //Post:        sets point3d location to (m_x, m_y, aZ)
  void setZ(const T aZ);
  
  //Description: Get z
  //Pre:         none
  //Post:        returns z location
  T z() const;
  
  point3d<T>& operator=(const point3d<T>& aRHS);
  template <class U>
  friend bool operator==(const point3d<U>& aLHS, const point3d<U>& aRHS);
  template <class U>
  friend bool operator!=(const point3d<U>& aLHS, const point3d<U>& aRHS);
  
  //Pre:         none
  //             class used in template needs to overload << operator
  //Post:        outputs aVector elements to aOutput
  //Description: output vector data
  template <class U>
  friend std::ostream& operator<<(std::ostream& aOutput, const point3d<U>& aPoint3d);
};

#include "point3d.hpp"

#endif /* defined(__hw6__point3d__) */
