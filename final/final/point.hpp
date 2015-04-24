//
//  Filename:     point.hpp
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  Here are the definitions of the functions declared in the point class
//


template <class T>
point<T>::~point()
{

}


template <class T>
void point<T>::set(const T aX, const T aY)
{
	setX(aX);
	setY(aY);
}


template <class T>
void point<T>::setX(const T aX)
{
	m_x = aX;
}


template <class T>
void point<T>::setY(const T aY)
{
	m_y = aY;
}


template <class T>
T point<T>::x() const
{
	return m_x;
}


template <class T>
T point<T>::y() const
{
	return m_y;
}
