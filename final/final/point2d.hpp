//
//  Filename:     point2d.hpp
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  Here are the definitions of the functions declared in the point2d class
//


template <class T>
point2d<T>::~point2d()
{

}


template <class T>
void point2d<T>::set(const T aX, const T aY)
{
	setX(aX);
	setY(aY);
}


template <class T>
void point2d<T>::setX(const T aX)
{
	m_x = aX;
}


template <class T>
void point2d<T>::setY(const T aY)
{
	m_y = aY;
}


template <class T>
T point2d<T>::x() const
{
	return m_x;
}


template <class T>
T point2d<T>::y() const
{
	return m_y;
}


template <class T>
point2d<T>& point2d<T>::operator=(const point2d<T>& aRHS)
{
	if (this != &aRHS)
	{
		setX(aRHS.x());
		setY(aRHS.y());
	}

	return *this;
}


template <class U>
bool operator==(const point2d<U>& aLHS, const point2d<U>& aRHS)
{
	return (aLHS.x() == aRHS.x() && aLHS.y() == aRHS.y());
}


template <class U>
bool operator!=(const point2d<U>& aLHS, const point2d<U>& aRHS)
{
	return !(aLHS == aRHS);
}


template <class U>
std::ostream& operator<<(std::ostream& aOutput, const point2d<U>& aPoint)
{
	std::stringstream ss;

	ss << "(" << aPoint.x() << ", " << aPoint.y() << ")";

	aOutput << ss.str();

	return aOutput;
}
