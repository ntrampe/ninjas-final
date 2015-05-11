//
//  Filename:     pde_notes.hpp
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  Here are the definitions of the functions declared in the
//                pde class
//


template <class T>
T pde_notes<T>::xLower(T aY) const
{
	return cos(M_PI*aY);
}


template <class T>
T pde_notes<T>::xUpper(T aY) const
{
	return exp(M_PI)*cos(M_PI*aY);
}


template <class T>
T pde_notes<T>::yLower(T aX) const
{
	return exp(M_PI*aX);
}


template <class T>
T pde_notes<T>::yUpper(T aX) const
{
	return -exp(M_PI*aX);
}

