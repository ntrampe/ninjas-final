//
//  Filename:     point3d.hpp
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  Here are the definitions of the functions declared in the point3d class
//


template <class T>
point3d<T>::~point3d()
{

}


template <class T>
void point3d<T>::set(const T aX, const T aY, const T aZ)
{
	this->set(aX, aY);
	setZ(aZ);
}


template <class T>
void point3d<T>::setZ(const T aZ)
{
	m_z = aZ;
}


template <class T>
T point3d<T>::z() const
{
	return m_z;
}


template <class T>
point3d<T>& point3d<T>::operator=(const point3d<T>& aRHS)
{
	if (this != &aRHS)
	{
		this->setX(aRHS.x());
		this->setY(aRHS.y());
		setZ(aRHS.z());
	}

	return *this;
}


template <class U>
bool operator==(const point3d<U>& aLHS, const point3d<U>& aRHS)
{
	return (aLHS.x() == aRHS.x() && aLHS.y() == aRHS.y() && aLHS.z() == aRHS.z());
}


template <class U>
bool operator!=(const point3d<U>& aLHS, const point3d<U>& aRHS)
{
	return !(aLHS == aRHS);
}


template <class U>
std::ostream& operator<<(std::ostream& aOutput, const point3d<U>& aPoint3d)
{
	std::stringstream ss;

	ss << "(" << aPoint3d.x() << ", " << aPoint3d.y() << ", " << aPoint3d.z() << ")";

	aOutput << ss.str();

	return aOutput;
}
