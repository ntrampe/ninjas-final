//
//  Filename:     pde_final.hpp
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  Here are the definitions of the functions declared in the
//                pde class
//


template <class T>
T pde_final<T>::xLower(T aY) const
{
	return sin(aY);
}


template <class T>
T pde_final<T>::xUpper(T aY) const
{
	return aY*0;
}


template <class T>
T pde_final<T>::yLower(T aX) const
{
	return sin(aX);
}


template <class T>
T pde_final<T>::yUpper(T aX) const
{
	return aX*0;
}

